#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2024 [SNALE - French SAS Company - RCS 951 724 616]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import inspect
import math
import os
from datetime import datetime
from datetime import timedelta

import cftime
import numpy as np
import pandas
from array_split import shape_split

from spatialetl.coverage.coverage import Coverage
from spatialetl.exception.not_found_in_rank_error import NotFoundInRankError
from spatialetl.operator.analytic.residence_time import compute_residence_time
from spatialetl.operator.interpolator.interpolator_core import temporal_2d_interpolation
from spatialetl.utils.logger import logging


class TimeCoverage(Coverage):
    """
    La classe TimeCoverage est une extension de la classe Coverage.
    Elle rajoute une dimension temporelle à la couverture horizontale classique.
    """

    TIME_DATUM = datetime(1970, 1, 1)
    TIME_DELTA = timedelta(minutes=15)
    TIME_INTERPOLATION_METHOD = "nearest"
    TIME_OVERLAPING_SIZE = 0

    def __init__(self, reader, bbox=None, resolution_x=None, resolution_y=None, start_time=None, end_time=None,
                 freq=None, nb_thread: int = os.cpu_count() - 1):

        Coverage.__init__(self, reader, bbox=bbox, resolution_x=resolution_x, resolution_y=resolution_y,
                          nb_thread=nb_thread);

        self.source_global_t_size = self.reader.get_t_size()
        self.source_global_axis_t = self.reader.read_axis_t(0, self.source_global_t_size, 0);

        self.temporal_resampling = False
        tmin = 0
        tmax = self.source_global_t_size
        zero_delta = timedelta(minutes=00)

        if start_time is not None:

            if type(start_time) == datetime:
                time = start_time
            elif type(start_time) == str:
                try:
                    time = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("start_time is not well formated " + str(ex))
            else:
                raise ValueError("start_time have to be string or datetime. Found " + str(type(start_time)))

            nearest_t_index = (np.abs(self.source_global_axis_t - time)).argmin()

            if time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                        '%Y-%m-%d %H:%M:%S') == zero_delta or abs(
                time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                         '%Y-%m-%d %H:%M:%S')) < TimeCoverage.TIME_DELTA:
                tmin = nearest_t_index
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                    TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        if end_time is not None:

            if type(end_time) == datetime:
                time = end_time
            elif type(end_time) == str:
                try:
                    time = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("end_time is not well formated " + str(ex))
            else:
                raise ValueError("end_time have to be string or datetime. Found " + str(type(end_time)))

            nearest_t_index = (np.abs(np.asarray(self.source_global_axis_t) - time)).argmin()

            if time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                        '%Y-%m-%d %H:%M:%S') == zero_delta or abs(
                time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                         '%Y-%m-%d %H:%M:%S')) < TimeCoverage.TIME_DELTA:
                tmax = nearest_t_index + 1
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                    TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        if freq is not None:
            self.temporal_resampling = True

            start_date = datetime.utcfromtimestamp(
                self.read_axis_t(type="source_global", with_overlap=False, timestamp=1)[tmin])
            end_date = datetime.utcfromtimestamp(
                self.read_axis_t(type="source_global", with_overlap=False,
                                 timestamp=1)[tmax - 1])

            self.target_global_axis_t = pandas.date_range(
                start=start_date,
                end=end_date,
                freq=freq).to_pydatetime()

            self.target_global_t_size = np.shape(self.target_global_axis_t)[0]

            if self.target_global_t_size == 0:
                raise ValueError("Target time axis is empty. Change start_time or end_time")

        else:
            self.target_global_axis_t = self.source_global_axis_t[tmin:tmax]
            self.target_global_t_size = tmax - tmin

        self.__init_parallel_map()

        if self.temporal_resampling and self.rank == 0:
            logging.info(
                '[temporal_interpolation] Source grid size : ' + str(self.source_global_t_size) + " time(s)")
            logging.info(
                '[temporal_interpolation] Target grid size : ' + str(self.target_global_t_size) + " time(s)")

        if self.horizontal_resampling:
            if self.rank == 0:
                logging.info(
                    '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                        self.source_global_y_size) + ")")
                logging.info(
                    '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                        self.target_global_y_size) + ")")

            self.compute_weight()

        if self.rank == 0 and self.comm:
            logging.debug("MPI map:")

        if self.comm:
            logging.debug(f"{'-' * 10} MPI rank n° {self.rank} {'-' * 10}")
        else:
            logging.debug("Multithreads map:")
            logging.debug(f"{'-' * 10} Source grid {'-' * 10}")
            for key in ['src_global_x', 'src_global_y', 'src_global_x_size', 'src_global_y_size', ]:
                logging.debug(f"{key} = {self.parallel_map[self.rank][key]}")
            logging.debug(f"{'-' * 10} Target grid {'-' * 10}")
            for key in ['dst_global_x', 'dst_global_y', 'dst_local_x_size', 'dst_local_y_size']:
                logging.debug(f"{key} = {self.parallel_map[self.rank][key]}")

        for key in self.parallel_map[self.rank]:
            if self.comm and key != "threads":
                logging.debug(f"{key} = {self.parallel_map[self.rank][key]}")

            # if len(self.parallel_map[self.rank]["threads"]) != 1:
        for thread in range(len(self.parallel_map[self.rank]["threads"])):
            logging.debug(f"   {'-' * 10} Rank {self.rank} - Thread n° {thread} {'-' * 10}")
            for thread_key in self.parallel_map[self.rank]["threads"][thread]:
                logging.debug(
                    f"    {thread_key} = {self.parallel_map[self.rank]['threads'][thread][thread_key]}")

        if self.rank == 0:
            logging.debug("-" * 20)

    def __init_parallel_map(self):
        """
        Create the MPI map for parallel processing.
        The MPI map is a dictionary that contains the mapping of the source and destination grids for each MPI rank.

        Examples
        --------
        >>> coverage.__init_parellel_map()
        """

        # Découpage des axes
        # if self.horizontal_resampling:
        #     # Découpage sur le temps uniquemenent
        #     target_slices = shape_split(target_sample, self.size, axis=[0, 1, 1])
        #     # Si on feet pas le nombre de proc
        #     if len( target_slices.flatten()) != self.size:
        #         target_slices = shape_split(target_sample, self.size, axis=[0, 0, 0])
        # else:

        target_mpi_sample = (self.target_global_t_size, self.target_global_y_size, self.target_global_x_size)

        # Split the axes with the MPI size
        target_mpi_slices = shape_split(target_mpi_sample, self.size, axis=[0, 0, 0])

        mpi_slice_index = 0
        for slyce in target_mpi_slices.flatten():
            mpi_slice = tuple(slyce)
            self.parallel_map[mpi_slice_index] = self.compute_slice_coordinates(
                mpi_slice,
                self.source_global_t_size,
                self.source_global_x_size,
                self.source_global_y_size,
                self.target_global_t_size,
                self.target_global_x_size,
                self.target_global_y_size
            )

            # Split the axes with the number of threads
            # We only split coverage grid
            target_threads_sample = (self.parallel_map[mpi_slice_index]["dst_local_y_size"],
                                     self.parallel_map[mpi_slice_index]["dst_local_x_size"])
            target_threads_slices = shape_split(target_threads_sample, self.threads_number, axis=[0, 0])

            # If we can't divide the grid dimensions with the number of threads,
            # we set the threads number with the number of slices
            self.threads_number = len(target_threads_slices.flatten())
            logging.debug(f"Threads number is {self.threads_number}")

            self.parallel_map[mpi_slice_index]['threads'] = np.empty([self.threads_number], dtype=object)

            slice_thread_index = 0
            for thread_slyce in target_threads_slices.flatten():
                thread_slice = tuple(thread_slyce)
                self.parallel_map[mpi_slice_index]["threads"][slice_thread_index] = Coverage.compute_slice_coordinates(
                    self,
                    thread_slice,
                    self.parallel_map[mpi_slice_index]["src_local_x_size"],
                    self.parallel_map[mpi_slice_index]["src_local_y_size"],
                    self.parallel_map[mpi_slice_index]["dst_local_x_size"],
                    self.parallel_map[mpi_slice_index]["dst_local_y_size"],
                    self.parallel_map[mpi_slice_index])
                slice_thread_index = slice_thread_index + 1

            mpi_slice_index = mpi_slice_index + 1

    def compute_slice_coordinates(self,
                                  slice,
                                  source_global_t_size: int,
                                  source_global_x_size: int,
                                  source_global_y_size: int,
                                  target_global_t_size: int,
                                  target_global_x_size: int,
                                  target_global_y_size: int,
                                  parent_slice=None
                                  ):
        """
        Compute slice coordinates in the source grid and the destination grid with overlap.

        Args:
           slice (array slice): Current slice to compute
           source_global_x_size (int) : Size of the global source x axis
           source_global_y_size (int) : Size of the global source y axis
           target_global_x_size (int) : Size of the global target x axis
           target_global_y_size (int) : Size of the global target y axis
           parent_slice (map of slice, optionial) : Slice of the parent slice
       """

        map = Coverage.compute_slice_coordinates(self,
                                                 slice[1:],
                                                 source_global_x_size,
                                                 source_global_y_size,
                                                 target_global_x_size,
                                                 target_global_y_size,
                                                 parent_slice
                                                 )

        # Grille source
        map["dst_global_t"] = slice[0]

        map["dst_local_t_size"] = map["dst_global_t"].stop - map["dst_global_t"].start

        dst_global_t_min_overlap = max(0, map["dst_global_t"].start - TimeCoverage.TIME_OVERLAPING_SIZE)
        dst_global_t_max_overlap = min(target_global_t_size,
                                       map["dst_global_t"].stop + TimeCoverage.TIME_OVERLAPING_SIZE)
        map["dst_global_t_overlap"] = np.s_[dst_global_t_min_overlap:dst_global_t_max_overlap]

        map["dst_global_t_size_overlap"] = map["dst_global_t_overlap"].stop - map["dst_global_t_overlap"].start

        dst_t_min = TimeCoverage.TIME_OVERLAPING_SIZE
        dst_t_max = map["dst_global_t_size_overlap"] - TimeCoverage.TIME_OVERLAPING_SIZE

        if map["dst_global_t"].start == 0:
            dst_t_min = 0

        if map["dst_global_t"].stop == target_global_t_size:
            dst_t_max = map["dst_global_t_size_overlap"]

        map["dst_local_t"] = np.s_[dst_t_min:dst_t_max]

        # Source grille
        if parent_slice is not None:
            source_global_axis_t = self.source_global_axis_t[parent_slice["dst_global_t"]]
        else:
            source_global_axis_t = self.source_global_axis_t

        target_global_axis_t = self.target_global_axis_t[map["dst_global_t"]]

        # TODO Add delta time
        idx = np.where(
            (source_global_axis_t >= np.min(
                target_global_axis_t - TimeCoverage.TIME_DELTA)) &
            (source_global_axis_t <= np.max(
                target_global_axis_t + TimeCoverage.TIME_DELTA)))

        tmin = np.min(idx[0])
        tmax = np.max(idx[0]) + 1

        # SRC GLOBAL
        map["src_global_t"] = np.s_[int(tmin):int(tmax)]
        map["src_global_t_size"] = tmax - tmin

        dst_global_t_min_overlap = max(0, map[
            "src_global_t"].start - TimeCoverage.TIME_OVERLAPING_SIZE)
        dst_global_t_max_overlap = min(source_global_t_size,
                                       map[
                                           "src_global_t"].stop + TimeCoverage.TIME_OVERLAPING_SIZE)
        map["src_global_t_overlap"] = np.s_[
            dst_global_t_min_overlap:dst_global_t_max_overlap]

        map["src_global_t_size_overlap"] = \
            map[
                "src_global_t_overlap"].stop - \
            map[
                "src_global_t_overlap"].start

        map["src_local_t_size"] = tmax - tmin
        map["src_local_t"] = np.s_[
            0:int(map["src_local_t_size"])]

        # OVERLAP
        map["src_local_t_size_overlap"] = \
            map[
                "src_global_t_size_overlap"]

        map["src_local_t_overlap"] = np.s_[
            0:map["src_local_t_size_overlap"]]

        return map

    # Axis
    def find_time_index(self, t, domain="source"):
        """Retourne l'index de la date la plus proche à TIME_DELTA_MIN prêt.
    @type t: datetime ou int
    @param t: date souhaitée ou l'index de la date souhaitée
    @return:  l'index de la date la plus proche à TIME_DELTA_MIN prêt ou une erreur si aucune date n'a pu être trouvée."""

        indexes_t = []

        if type(t) == int or type(t) == np.int32 or type(t) == np.int64:

            if t < 0 or t >= self.get_t_size():
                raise ValueError("Time index have to range between 0 and " + str(
                    self.get_t_size() - 1) + ". Actually Time index = " + str(t))

            indexes_t.append(int(t));

        elif type(t) == datetime or type(t) == cftime._cftime.datetime or type(t) == cftime._cftime.real_datetime or type(t) == cftime._cftime.DatetimeGregorian:

            target_timestamp = (t - TimeCoverage.TIME_DATUM).total_seconds()
            array = self.read_axis_t(type="source_mpi", timestamp=1)

            logging.debug("[TimeCoverage][find_time_index()] Looking for : " + str(t))

            X = np.abs(array - target_timestamp)
            idx = np.where(X <= (TimeCoverage.TIME_DELTA).total_seconds())

            if (len(idx[0]) == 1):
                index_t = idx[0][0]
                if domain == "source":
                    indexes_t.append(int(index_t))
                elif domain == "source_global":
                    indexes_t.append(self.parallel_map[self.rank]["src_global_t"].start + int(index_t))
                else:
                    raise ValueError("Type doesn't match [source, source_global]")

                logging.debug(
                    f"[TimeCoverage][find_time_index()] Found : {self.read_axis_t(type='source_mpi', timestamp=0)[int(index_t)]}")

            else:
                for index in range(np.shape(idx)[1]):
                    index_t = idx[0][index]

                    if domain == "source":
                        indexes_t.append(int(index_t))
                    elif domain == "source_global":
                        indexes_t.append(self.parallel_map[self.rank]["src_global_t"].start + int(index_t))
                    else:
                        raise ValueError("Type doesn't match [source, source_global]")

                    logging.debug(
                        f"[TimeCoverage][find_time_index()] Found : {self.read_axis_t(type='source_mpi', timestamp=0)[int(index_t)]}")

            if not indexes_t:
                raise NotFoundInRankError(self.rank,
                                          "'" + str(
                                              t) + "' not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                                              TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        else:
            raise ValueError("" + str(t) + " have to be an integer or a datetime. Current type: " + str(type(t)))

        indexes_t = np.array(np.sort(np.unique(indexes_t)))
        logging.debug("[TimeCoverage][find_time_index()] Found " + str(len(indexes_t)) + " candidate datetime(s)")

        return indexes_t

    def read_axis_t(self, type="target_mpi", with_overlap=False, timestamp=0, current_thread: int = None):
        """Retourne les valeurs de l'axe t.
    @param timestamp: égale 1 si le temps est souhaité en timestamp depuis TIME_DATUM.
    @return:  un tableau à une dimensions [z] au format datetime ou timestamp si timestamp=1."""
        if type == "target_global":
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t];
            return self.target_global_axis_t

        elif type == "source_global":
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.source_global_axis_t];
            return self.source_global_axis_t

        elif type == "source_mpi" and with_overlap is True:
            return self.reader.read_axis_t(self.parallel_map[self.rank]["src_global_t_overlap"].start,
                                           self.parallel_map[self.rank]["src_global_t_overlap"].stop, timestamp)

        elif type == "source_mpi" and with_overlap is False:
            return self.reader.read_axis_t(self.parallel_map[self.rank]["src_global_t"].start,
                                           self.parallel_map[self.rank]["src_global_t"].stop, timestamp)

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_t(
                self.parallel_map[self.rank]["threads"][current_thread]["src_global_t_overlap"].start,
                self.parallel_map[self.rank]["threads"][current_thread]["src_global_t_overlap"].stop, timestamp)

        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_t(
                self.parallel_map[self.rank]["threads"][current_thread]["src_global_t"].start,
                self.parallel_map[self.rank]["threads"][current_thread]["src_global_t"].stop, timestamp)

        elif type == "target" and with_overlap is True:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[
                            self.parallel_map[self.rank]["threads"][current_thread]["dst_global_t_overlap"]]];
            return self.target_global_axis_t[
                self.parallel_map[self.rank]["threads"][current_thread]["dst_global_t_overlap"]]

        elif type == "target" and with_overlap is False:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[
                            self.parallel_map[self.rank]["threads"][current_thread]["dst_global_t"]]];
            return self.target_global_axis_t[self.parallel_map[self.rank]["threads"][current_thread]["dst_global_t"]]

        elif type == "target_mpi" and with_overlap is True:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t_overlap"]]];
            return self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t_overlap"]]

        else:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t"]]];
            return self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t"]]

    def get_t_size(self, type="target", with_overlap=False):
        if type == "target_global":
            return self.target_global_t_size
        elif type == "source_global":
            return self.source_global_t_size
        elif type == "source":
            return self.parallel_map[self.rank]["src_local_t_size"]
        elif type == "target" and with_overlap is True:
            return self.parallel_map[self.rank]["dst_local_t_size_overlap"]
        else:
            return self.parallel_map[self.rank]["dst_local_t_size"]

    def __read_variable(self, function_name, time, is_binary_data=False):

        fn = getattr(self.reader, function_name)

        index_t = self.find_time_index(time);

        layers = np.stack([fn(
            self.parallel_map[self.rank]["src_global_t"].start + index_t[t],
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop) for t in range(0, len(index_t))])

        times = self.read_axis_t(type="source_mpi")
        source_time = times[index_t]

        data = temporal_2d_interpolation(source_time, time, layers, TimeCoverage.TIME_INTERPOLATION_METHOD,
                                         TimeCoverage.TIME_DELTA)

        is_vector = True if len(np.shape(data)) == 3 and np.shape(data)[0] == 2 else False

        if self.horizontal_resampling:
            if is_vector:
                return [
                    self.resample_2d_variable(data[0], is_binary_data=is_binary_data),
                    self.resample_2d_variable(data[1], is_binary_data=is_binary_data)
                ]
            else:
                return self.resample_2d_variable(data, is_binary_data=is_binary_data)
        else:
            if is_vector:
                return [
                    data[0][self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]],
                    data[1][self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]
                ]
            else:
                return data[self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]

    # Variables
    def read_variable_2D_sea_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t, is_binary_data=True)

    def read_variable_2D_wet_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t, is_binary_data=True)

    def read_variable_2D_land_binary_mask_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t, is_binary_data=True)

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_height_above_geoid_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_column_thickness_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_temperature_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_salinity_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_pressure_at_time(self, t):
        """Retourne la pression à la surface de la mer (sea surface pressure) à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_density_at_time(self, t):
        """Retourne la densité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_turbidity_at_time(self, t):
        """Retourne la turbidité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_residence_time_at_time(self, t):
        """
            Read the sea surface residence time at a given time
            @param t: date souhaitée
            @return: un tableau en deux dimensions [y,x].
        """
        fn = getattr(self.reader, "read_variable_sea_surface_biogeochemical_tracer_1_at_time")

        indexes_t = self.find_time_index(t);

        layers = np.stack([fn(
            self.parallel_map[self.rank]["src_global_t"].start + indexes_t[index_t],
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop) for index_t in range(0, len(indexes_t))])

        times = self.read_axis_t(type="source_mpi")
        source_time = times[indexes_t] - np.min(times)

        def to_decimal_days(td):
            return td.total_seconds() / (3600 * 24)

        vectorized_square = np.vectorize(to_decimal_days)
        result = vectorized_square(source_time)

        data = compute_residence_time(result, layers)

        if self.horizontal_resampling:
            return self.resample_2d_variable(data)
        else:
            return data[self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]

    #################
    # HYDRO
    # Ground level
    #################

    def read_variable_sea_water_temperature_at_ground_level_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_salinity_at_ground_level_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_velocity_at_ground_level_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
           @type t: datetime ou l'index
           @param t: date souhaitée
           @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_residence_time_at_ground_level_at_time(self, t):
        """
            Read the residence time at ground level at a given time
            @param t: date souhaitée
            @return: un tableau en deux dimensions [y,x].
        """
        fn = getattr(self.reader, "read_variable_biogeochemical_tracer_1_at_ground_level_at_time")

        indexes_t = self.find_time_index(t);

        layers = np.stack([fn(
            self.parallel_map[self.rank]["src_global_t"].start + indexes_t[index_t],
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop) for index_t in range(0, len(indexes_t))])

        times = self.read_axis_t(type="source_mpi")
        source_time = times[indexes_t] - np.min(times)

        def to_decimal_days(td):
            return td.total_seconds() / (3600 * 24)

        vectorized_square = np.vectorize(to_decimal_days)
        result = vectorized_square(source_time)

        data = compute_residence_time(result, layers)

        if self.horizontal_resampling:
            return self.resample_2d_variable(data)
        else:
            return data[self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_barotropic_sea_water_velocity_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_barotropic_sea_water_speed_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)
        return np.hypot(comp[0], comp[1])

    def read_variable_barotropic_sea_water_from_direction_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) + 180.0 % 360.0

        return result

    def read_variable_barotropic_sea_water_to_direction_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) % 360.0

        return result

    def read_variable_ocean_tracer_residence_time_at_time(self, t):
        """ Read ocean tracer residence time at the given time
           @type t: datetime ou l'index
           @param t: date souhaitée
           @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self, t):
        """Retourne la hauteur significative des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_breaking_height_at_time(self, t):
        """Retourne la hauteur de déferlement des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_mean_period_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_peak_period_at_time(self, t):

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_from_direction_at_time(self, t):

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_to_direction_at_time(self, t):

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self, t):
        """Retourne la dérive de Stokes en surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_radiation_pressure_bernouilli_head_at_time(self, t):
        """Retourne la pression J due aux vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self, t):
        """Retourne la waves_to_ocean_energy_flux à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self, t):
        """Retourne la l'énergie des vagues dissipée par le fond à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    #################
    # WAVES
    # Momentum flux
    #################
    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self, t):
        """Retourne la composante u du tau atmosphere->vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_waves_momentum_flux_to_ocean_at_time(self, t):
        """Retourne la composante u du tau vagues->ocean à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    #################
    # METEO
    # 2D
    #################
    def read_variable_rainfall_amount_at_time(self, t):
        """Retourne les composantes u,v de rain à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_air_pressure_at_time(self, t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_wind_stress_at_time(self, t):
        """Retourne les composantes u,v de la contrainte de vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface sensible heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_downward_latent_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface latente heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_air_temperature_at_time(self, t):
        """Retourne les composantes u,v de surface air temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_dew_point_temperature_at_time(self, t):
        """Retourne les composantes u,v de dewpoint temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_downward_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_downward_thermal_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_thermal_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self, t):
        """Retourne les composantes u,v du vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_wind_speed_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        return np.hypot(comp[0], comp[1])

    def read_variable_wind_from_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) + 180.0 % 360.0

        return result

    def read_variable_wind_to_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) % 360.0

        return result
