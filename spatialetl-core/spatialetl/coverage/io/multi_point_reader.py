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
import numpy as np
from scipy.interpolate import griddata, NearestNDInterpolator

from spatialetl.coverage import TimeCoverage
from spatialetl.coverage.io.memory_reader import MemoryReader
from spatialetl.operator.interpolator.smoothing import nanmoving_average2
from spatialetl.point import TimeMultiPoint
from spatialetl.utils.logger import logging


class MultiPointReader(MemoryReader):

    def __init__(self,
                 points,
                 resolution_x,
                 resolution_y,
                 mask_axis_x=None,
                 mask_axis_y=None,
                 mask_values=None
                 ):

        self.points = points
        x = points.read_axis_x()
        y = points.read_axis_y()
        self.target_axis_x = np.arange(min(x), max(x), resolution_x)
        self.target_axis_y = np.arange(min(y), max(y), resolution_y)
        self.target_grid_x, self.target_grid_y = np.meshgrid(self.target_axis_x, self.target_axis_y)
        self.src_points = (points.read_axis_x().flatten(), points.read_axis_y().flatten())

        if mask_axis_x is not None and mask_axis_y is not None and mask_values is not None:
            # If needed, we interpole the land sea mask on the destination grid
            logging.info("[MultiPointReader] Interpolate sea binary mask")
            src_grid_x, src_grid_y = np.meshgrid(mask_axis_x, mask_axis_y)
            points = np.array([src_grid_x.flatten(), src_grid_y.flatten()]).T
            interp = NearestNDInterpolator(points, mask_values.flatten())
            self.mask = interp((self.target_grid_x, self.target_grid_y))
        else:
            self.mask = None

        if type(self.points) == TimeMultiPoint:
            self.t = self.points.read_axis_t()

    def is_regular_grid(self):
        return True

    def get_t_size(self):
        return np.shape(self.t)[0];

    def get_x_size(self):
        return np.shape(self.target_axis_x)[0];

    def get_y_size(self):
        return np.shape(self.target_axis_y)[0];

    def read_axis_x(self, xmin, xmax, ymin=None, ymax=None):
        return self.target_axis_x[xmin:xmax]

    def read_axis_y(self, xmin=None, xmax=None, ymin=None, ymax=None):
        return self.target_axis_y[ymin:ymax]

    def read_axis_t(self, tmin, tmax, timestamp):
        if timestamp == 1:
            return np.asarray([(t - TimeCoverage.TIME_DATUM).total_seconds() \
                               for t in self.t[tmin:tmax]]);
        else:
            return self.t[tmin:tmax]

    # Variables
    def read_variable_longitude(self, xmin, xmax, ymin, ymax):
        return self.read_axis_x(xmin, xmax, ymin, ymax)

    def read_variable_latitude(self, xmin, xmax, ymin, ymax):
        return self.read_axis_y(xmin, xmax, ymin, ymax)

    def read_variable_2D_land_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_2D_sea_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_3D_sea_binary_mask_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_3D_land_binary_mask_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    #################
    # HYDRO
    # 2D
    #################

    def read_variable_ocean_tracer_residence_time_at_time(self, index_t, xmin, xmax, ymin, ymax):
        # We get the datetime of the index_t because we need to find all times records available according to TimeMultiPoint.TIME_DELTA
        time = self.t[index_t]
        values = self.points.read_variable_ocean_tracer_residence_time_at_time(time)
        data = griddata(self.src_points, values.flatten(), (self.target_grid_x, self.target_grid_y), method='linear')

        if self.mask is not None:
            data[self.mask != 1] = np.nan

        data = nanmoving_average2(data, 6)
        return data[ymin:ymax, xmin:xmax]

    #################
    # METEO
    # 2D
    #################

    def read_variable_rainfall_amount_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, index_t, xmin, xmax, ymin, ymax):
        return self.sp[index_t, ymin:ymax, xmin:xmax]

    def read_variable_sea_surface_air_pressure_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_air_temperature_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_dew_point_temperature_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_solar_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_thermal_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_solar_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_thermal_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")
