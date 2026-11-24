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

import numpy as np
from datetime import datetime

from spatialetl.exception.variable_name_error import VariableNameError
from spatialetl.providers.symphonie.coverage.netcdf.symphonie_reader import SYMPHONIEReader as AbstractSYMPHONIEReader
from spatialetl.utils.logger import logging


class SYMPHONIEEulerianTracersReader(AbstractSYMPHONIEReader):
    """
La classe SymphonieReader permet de lire les données du format Symphonie

@param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
@param myFile : lien vers le fichier de données (que l'on trouve dans GRAPHIQUES)
"""

    def __init__(self, myGrid, myFile, tracer_index=1):
        AbstractSYMPHONIEReader.__init__(self, myGrid, myFile);

        self.tracer_index = tracer_index

    def set_tracer_index(self, index):
        self.tracer_index = index

    def find_first_nonzero(self, variable, index_z, index_start, index_end):
        """
           Find the first time index at which a tracer contains non-zero values.

           A binary search is used to reduce the number of NetCDF files that
           need to be opened and inspected. The search assumes that once a
           tracer becomes non-zero, it remains non-zero at subsequent time steps.

           Args:
               variable (str): Name of the NetCDF variable to inspect (e.g. 'bio1').
               index_z (int): Vertical index at which the variable is evaluated.
               index_start (int): First time index of the search interval.
               index_end (int): Last time index of the search interval.

           Returns:
               int: First time index at which at least one finite, non-zero value
                   is found. If no non-zero value is found in the search interval,
                   the returned index corresponds to the last index checked.

           Note:
               The binary search assumes a monotonic transition from zero to
               non-zero values. It may return incorrect results if the variable
               becomes zero again after its first non-zero occurrence.
           """

        while index_start < index_end:

            index_mid = (index_start + index_end) // 2

            self.open_file(index_mid)

            data = self.ncfile.variables[variable][0, index_z, :, :]

            if np.any((data != 0) & np.isfinite(data)):
                index_end = index_mid
            else:
                index_start = index_mid + 1

        return index_start

    def get_date_tracers_map(self):
        """
           Build a dictionary associating each tracer with its first appearance date.

           Tracers are identified by sequential NetCDF variables named 'bio1',
           'bio2', 'bio3', etc. The search stops when the next tracer variable
           is not found.

           For each tracer, the first time index containing at least one finite,
           non-zero value is determined using binary search through
           `find_first_nonzero()`.

           Since tracers are assumed to appear sequentially in time, the search
           for each tracer starts immediately after the first appearance of
           the previous tracer, reducing the number of files to inspect.

           Returns:
               dict[int, datetime]: Dictionary mapping each tracer number to its
                   first appearance date.

                   If no non-zero value is found for a tracer, the first date
                   in `self.times` is used as a fallback.

           Note:
               - The NetCDF files must be ordered chronologically.
               - Tracer appearance dates must be strictly increasing.
               - Once a tracer becomes non-zero, it must remain non-zero at
                 subsequent time steps.
               - Dates are converted from cftime.datetime to datetime.datetime.
           """
        date_tracers_map = {}
        index_z = self.get_z_size() - 1

        index_start = 0
        index_end = len(self.times) - 1

        n = 1
        while f"bio{n}" in self.ncfile.variables:

            variable = f"bio{n}"

            index_found = self.find_first_nonzero(
                variable,
                index_z,
                index_start,
                index_end
            )

            self.open_file(index_found)

            data = self.ncfile.variables[variable][0, index_z, :, :]

            if np.any((data != 0) & np.isfinite(data)):
                date_tracers_map[n] = datetime(
                    self.times[index_found].year, self.times[index_found].month, self.times[index_found].day,
                    self.times[index_found].hour, self.times[index_found].minute, self.times[index_found].second,
                    self.times[index_found].microsecond
                )

                index_start = index_found + 1
            else:
                date_tracers_map[n] = datetime(
                    self.times[0].year, self.times[0].month, self.times[0].day,
                    self.times[0].hour, self.times[0].minute, self.times[0].second,
                    self.times[0].microsecond
                )

            n += 1

        return date_tracers_map

    #################
    # HYDRO
    # Sea Surface
    #################

    def read_variable_sea_surface_biogeochemical_tracer_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            if f"bio{self.tracer_index}" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables[f"bio{self.tracer_index}"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug(f"No variables found for 'sea_surface_biogeochemical_tracer_{self.tracer_index}'")
                raise (VariableNameError("SymphonieEulerianTracersReader",
                                         f"No variables found for 'sea_surface_biogeochemical_tracer_{self.tracer_index}'",
                                         1000))

            if AbstractSYMPHONIEReader.APPLY_WET_MASK and "wetmask_t" in self.ncfile.variables:  # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieEulerianTracersReader", "An error occured : '" + str(ex) + "'", 1000))

        #################
        # HYDRO
        # Ground level
        #################

    def read_variable_biogeochemical_tracer_at_ground_level_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = 0
            if f"bio{self.tracer_index}" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables[f"bio{self.tracer_index}"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug(f"No variables found for 'biogeochemical_tracer_{self.tracer_index}_at_ground_level'")
                raise (VariableNameError("SymphonieEulerianTracersReader",
                                         f"No variables found for 'biogeochemical_tracer_{self.tracer_index}_at_ground_level'",
                                         1000))

            if AbstractSYMPHONIEReader.APPLY_WET_MASK and "wetmask_t" in self.ncfile.variables:  # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieEulerianTracersReader", "An error occured : '" + str(ex) + "'", 1000))
