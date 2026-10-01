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
from __future__ import division, print_function, absolute_import

import glob
import os
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from netCDF4 import Dataset, num2date

from spatialetl.exception.variable_name_error import VariableNameError
from spatialetl.point import TimeMultiPoint
from spatialetl.point.io.multi_point_reader import MultiPointReader
from spatialetl.providers.symphonie.coverage.netcdf.symphonie_reader import extract_times_from_file, normalize_units
from spatialetl.utils.logger import logging
from spatialetl.utils.variable_definition import VariableDefinition


def extract_times_from_file(file):
    times = []

    try:
        with Dataset(file, 'r') as ds:
            time_var = ds.variables.get("time")
            if time_var is None:
                return times

            units = normalize_units(getattr(time_var, "long_name").replace('current time in', ''))
            calendar = getattr(time_var, "calendar", "gregorian")

            current_time = time_var[:]

            # vectorized conversion when possible
            if current_time.shape == (1,):
                times.append(num2date(current_time[0], units=units, calendar=calendar))
            else:
                decoded = num2date(current_time, units=units, calendar=calendar)
                times.extend(t.replace(microsecond=0) for t in decoded)

    except Exception as ex:
        raise ValueError(f"Unable to decode time records in file {file}: {ex}")
    finally:
        return times


class SYMPHONIELagrangianDriftersReader(MultiPointReader):
    def __init__(self, myFile):
        """
        """
        initial_drifters_filename = os.path.join(myFile, "drifter_initial.nc")
        MultiPointReader.__init__(self, initial_drifters_filename)

        # Read initial drifters
        try:
            initial_nc = Dataset(initial_drifters_filename, 'r')
            self.names = np.ma.filled(initial_nc.variables["init_identifier"][:], np.nan)
            self.x = np.ma.filled(initial_nc.variables["init_lon"][:], np.nan)
            self.y = np.ma.filled(initial_nc.variables["init_lat"][:], np.nan)
            self.z = np.ma.filled(initial_nc.variables["init_z"][:], np.nan)
            self.ref_identifiers = self.names.astype(int)
        except Exception as e:
            logging.error(e)
            logging.debug(e)

        if myFile is not None:
            if os.path.isfile(myFile):
                self.files = [myFile]
            elif os.path.isdir(myFile):
                self.files = sorted(glob.glob(os.path.join(myFile, "*.nc")))
            elif myFile.endswith("*"):
                self.files = sorted(glob.glob(myFile + ".nc"))
            else:
                raise ValueError("Unable to decode file " + str(myFile))

            self.ncfile = Dataset(self.files[0], 'r')
            self.last_opened_t_index = 0
        else:
            self.files = []
            self.ncfile = Dataset(self.filename, 'r')
            self.last_opened_t_index = 0

        self.t_size = len(self.files)
        self.times = []

        if self.t_size > 20:
            # Multithreading process
            with ProcessPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
                futures = [executor.submit(extract_times_from_file, f) for f in self.files]

                for future in as_completed(futures):
                    result = future.result()
                    self.times.extend(result)
        else:
            # Sequential process
            for f in self.files:
                results = extract_times_from_file(f)
                self.times.extend(results)

        # Note : sort in times and files has to be exactly the same order
        self.t_size = len(self.times)

        if len(self.times) == 0:
            logging.info("No time records found")

    def open_file(self, index_t):

        if len(self.files) == 1 and index_t > 0:
            return

        if index_t != self.last_opened_t_index:
            self.close()
            self.ncfile = Dataset(self.files[index_t])
            self.last_opened_t_index = index_t

    def close(self):
        self.ncfile.close()

    # Axis
    def get_x_size(self):
        return len(self.x);

    def get_y_size(self):
        return len(self.y);

    def get_z_size(self):
        return len(self.z);

    def get_t_size(self):
        return len(self.times);

    def read_axis_x(self):
        return self.x

    def read_axis_y(self):
        return self.y

    def read_axis_z(self):
        return self.z

    def read_axis_t(self, tmin, tmax, timestamp=0):
        """
        Returns the time axis
        :param tmin: Min slice index
        :param tmax: Max slice index
        :param timestamp: True if time in timestamp, else datime
        :return: Array of time
        """

        if timestamp == 1:
            return np.asarray([(t - TimeMultiPoint.TIME_DATUM).total_seconds() \
                               for t in self.times[tmin:tmax]]);
        else:
            return np.asarray(self.times[tmin:tmax])

    def read_metadata(self):
        metadata = {};

        return metadata

    def read_variable_point_names(self):
        return self.names

    def read_variable_ocean_xy_tracer_at_time(self, index_t):
        try:
            self.open_file(index_t)

            if "latitude" in self.ncfile.variables and "longitude" in self.ncfile.variables and "identifier" in self.ncfile.variables:
                current_identifiers = np.ma.filled(self.ncfile.variables["identifier"][0], fill_value=-9999).astype(int)

                # Remove values for drifters outside the domain
                mask = current_identifiers <= 0.
                lat = np.ma.filled(self.ncfile.variables["latitude"][0], fill_value=np.nan)
                lon = np.ma.filled(self.ncfile.variables["longitude"][0], fill_value=np.nan)
                lat[mask] = np.nan
                lon[mask] = np.nan

                # Re-order values with the initial identifiers order
                positions = defaultdict(deque)
                for i, val in enumerate(current_identifiers):
                    positions[val].append(i)

                indices = np.fromiter(
                    (positions[val].popleft() if positions[val] else -1
                     for val in self.ref_identifiers),
                    dtype=int,
                    count=len(self.ref_identifiers)
                )

                # Extract data and make a 2D array (lat, lon)
                data = np.vstack((lat[indices], lon[indices]))

            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['ocean_tracer_residence_time']) + "'")
                raise (VariableNameError("SYMPHONIELagrangianDriftersReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME[
                                                 'ocean_tracer_residence_time']) + "'",
                                         1000))

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SYMPHONIELagrangianDriftersReader",
                                     "An error occured at time index '" + str(index_t) + "' : '" + str(ex) + "'", 1000))
