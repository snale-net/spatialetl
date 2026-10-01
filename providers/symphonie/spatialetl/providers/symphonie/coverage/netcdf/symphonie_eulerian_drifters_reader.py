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

from spatialetl.providers.symphonie.coverage.netcdf.symphonie_reader import SYMPHONIEReader as AbstractSYMPHONIEReader
from spatialetl.exception.variable_name_error import VariableNameError
from spatialetl.utils.variable_definition import VariableDefinition
from spatialetl.utils.logger import logging


class SYMPHONIEEulerianDriftersReader(AbstractSYMPHONIEReader):
    """
La classe SymphonieReader permet de lire les données du format Symphonie

@param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
@param myFile : lien vers le fichier de données (que l'on trouve dans GRAPHIQUES)
"""

    def __init__(self,myGrid, myFile):
        AbstractSYMPHONIEReader.__init__(self, myGrid, myFile);

    #################
    # HYDRO
    # Sea Surface
    #################

    def read_variable_sea_surface_biogeochemical_tracer_1_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            if "bio1" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["bio1"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'sea_surface_biogeochemical_tracer_1'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'sea_surface_biogeochemical_tracer_1'",
                                         1000))

            if AbstractSYMPHONIEReader.APPLY_WET_MASK and "wetmask_t" in self.ncfile.variables:  # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # HYDRO
    # Ground level
    #################

    def read_variable_biogeochemical_tracer_1_at_ground_level_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = 0
            if "bio1" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["bio1"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'biogeochemical_tracer_1_at_ground_level'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'biogeochemical_tracer_1_at_ground_level'",
                                         1000))

            if AbstractSYMPHONIEReader.APPLY_WET_MASK and "wetmask_t" in self.ncfile.variables:  # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))