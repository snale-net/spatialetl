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
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import int16, float32, float64

from spatialetl.coverage.io.coverage_writer import CoverageWriter
from spatialetl.coverage.level_coverage import LevelCoverage
from spatialetl.coverage.time_coverage import TimeCoverage
from spatialetl.coverage.time_level_coverage import TimeLevelCoverage
from spatialetl.exception.coverage_error import CoverageError
from spatialetl.utils.logger import logging
from spatialetl.utils.variable_definition import VariableDefinition


class DefaultWriter(CoverageWriter):

    def __init__(self, cov, myFile, mode='w'):
        CoverageWriter.__init__(self, cov, myFile);
        self.mode = mode
        format = 'NETCDF4_CLASSIC'

        if self.mode == 'w' and self.coverage.rank == 0:
            self.ncfile = Dataset(self.filename, 'w', format=format)
            self.ncfile.description = 'Generated with pySpatialETL'

            # dimensions
            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['latitude'],
                                        self.coverage.get_y_size(type="target_global"))
            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['longitude'],
                                        self.coverage.get_x_size(type="target_global"))

            if self.coverage.is_regular_grid() == True:

                # variables
                latitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['latitude'], float32,
                                                       (VariableDefinition.VARIABLE_NAME['latitude'],))
                latitudes.long_name = VariableDefinition.LONG_NAME['latitude'];
                latitudes.standard_name = VariableDefinition.STANDARD_NAME['latitude'];
                latitudes.valid_min = -90.;
                latitudes.valid_max = 90.;
                latitudes.axis = "Y";
                latitudes.units = VariableDefinition.CANONICAL_UNITS['latitude'];

                longitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['longitude'], float32,
                                                        (VariableDefinition.VARIABLE_NAME['longitude'],))
                longitudes.long_name = VariableDefinition.LONG_NAME['longitude'];
                longitudes.standard_name = VariableDefinition.STANDARD_NAME['longitude'];
                longitudes.valid_min = -180.;
                longitudes.valid_max = 180.;
                longitudes.axis = "X";
                longitudes.units = VariableDefinition.CANONICAL_UNITS['longitude'];

                latitudes[self.coverage.parallel_map[self.coverage.rank]["dst_global_y"]] = self.coverage.read_axis_y()
                longitudes[self.coverage.parallel_map[self.coverage.rank]["dst_global_x"]] = self.coverage.read_axis_x()

            else:

                latitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['latitude'], float32,
                                                       (VariableDefinition.VARIABLE_NAME['latitude'],
                                                        VariableDefinition.VARIABLE_NAME['longitude'],))
                latitudes.long_name = VariableDefinition.LONG_NAME['latitude'];
                latitudes.standard_name = VariableDefinition.STANDARD_NAME['latitude'];
                latitudes.valid_min = -90.;
                latitudes.valid_max = 90.;
                latitudes.axis = "Y";
                latitudes.units = VariableDefinition.CANONICAL_UNITS['latitude'];

                longitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['longitude'], float32,
                                                        (VariableDefinition.VARIABLE_NAME['latitude'],
                                                         VariableDefinition.VARIABLE_NAME['longitude'],))
                longitudes.long_name = VariableDefinition.LONG_NAME['longitude'];
                longitudes.standard_name = VariableDefinition.STANDARD_NAME['longitude'];
                longitudes.valid_min = -180.;
                longitudes.valid_max = 180.;
                longitudes.axis = "X";
                longitudes.units = VariableDefinition.CANONICAL_UNITS['longitude'];

                # data
                latitudes[
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_y"],
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_x"]] = self.coverage.read_axis_y()

                longitudes[
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_y"],
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_x"]] = self.coverage.read_axis_x()

            if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):
                self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['time'],
                                            self.coverage.get_t_size(type="target_global"))
                times = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['time'], float64,
                                                   (VariableDefinition.VARIABLE_NAME['time'],))
                times.units = 'seconds since 1970-01-01 00:00:00'
                times.calendar = 'gregorian'
                times.standard_name = 'time'
                times.axis = 'T'
                times.conventions = "UTC time"

                times[
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_t"]] = date2num(
                    self.coverage.read_axis_t(), units=times.units, calendar=times.calendar)

            if (isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

                if self.coverage.is_sigma_coordinate():
                    raise ValueError(
                        "This writer supports only Coverage with a regular vertical axis. Use 'resolution_z' on the Coverage to interpolate")

                self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['depth'],
                                            self.coverage.get_z_size(type="target"))
                levels = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['depth'], float64,
                                                    (VariableDefinition.VARIABLE_NAME['depth'],))
                levels.standard_name = 'depth'
                levels.long_name = "Positive depth"
                levels.axis = 'Z'
                levels.units = VariableDefinition.CANONICAL_UNITS['depth'];

                levels[:] = self.coverage.read_axis_z()

        if self.mode == "a":

            self.ncfile = Dataset(self.filename, 'a', format=format)

            # dimensions
            if self.ncfile.dimensions['latitude'].size != self.coverage.get_y_size(type="target_global"):
                raise ValueError(
                    "Latitude dimension hasn't the same size than the Coverage. Unable to append the file.")

            if self.ncfile.dimensions['longitude'].size != self.coverage.get_x_size(type="target_global"):
                raise ValueError(
                    "Longitude dimension hasn't the same size than the Coverage. Unable to append the file.")

            if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

                if self.ncfile.dimensions['time'].size != self.coverage.get_t_size(type="target_global"):
                    raise ValueError(
                        "Time dimension hasn't the same size than the Coverage. Unable to append the file.")

            if (isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

                if self.coverage.is_sigma_coordinate():
                    raise ValueError("This writer supports only Coverage with a regular vertical axis.")

                if self.ncfile.dimensions['depth'].size != self.coverage.get_z_size(type="target"):
                    raise ValueError(
                        "Depth dimensions hasn't the same size than the Coverage. Unable to append the file.")

    def close(self):
        if self.coverage.rank == 0:
            self.ncfile.close()

    def _create_nc_var(self, variable_name,
                       variable_long_name,
                       variable_standard_name,
                       variable_units,
                       variable_comment=None,
                       is_time_var=False,
                       is_depth_var=False,
                       variable_type=float32,
                       variable_fill_value=9.96921e+36):
        if self.coverage.rank == 0:
            if variable_name in self.ncfile.variables:
                var = self.ncfile.variables[variable_name]
            else:
                if is_time_var and is_depth_var:
                    var = self.ncfile.createVariable(variable_name, variable_type,
                                                     (VariableDefinition.VARIABLE_NAME['time'],
                                                      VariableDefinition.VARIABLE_NAME['depth'],
                                                      VariableDefinition.VARIABLE_NAME['latitude'],
                                                      VariableDefinition.VARIABLE_NAME['longitude'],),
                                                     fill_value=variable_fill_value)
                elif is_depth_var:
                    var = self.ncfile.createVariable(variable_name, variable_type,
                                                     (VariableDefinition.VARIABLE_NAME['depth'],
                                                      VariableDefinition.VARIABLE_NAME['latitude'],
                                                      VariableDefinition.VARIABLE_NAME['longitude'],),
                                                     fill_value=variable_fill_value)
                elif is_time_var:
                    var = self.ncfile.createVariable(
                        variable_name,
                        variable_type,
                        (VariableDefinition.VARIABLE_NAME['time'],
                         VariableDefinition.VARIABLE_NAME['latitude'],
                         VariableDefinition.VARIABLE_NAME['longitude'],),
                        fill_value=variable_fill_value)
                else:
                    var = self.ncfile.createVariable(variable_name, variable_type,
                                                     (VariableDefinition.VARIABLE_NAME['latitude'],
                                                      VariableDefinition.VARIABLE_NAME['longitude'],),
                                                     fill_value=variable_fill_value)

            var.long_name = variable_long_name
            var.standard_name = variable_standard_name
            var.units = variable_units

            if variable_comment:
                var.comment = variable_comment
        else:
            var = None

        return var

    def _write_2d_scalar(self,
                         function_name,
                         variable_name,
                         variable_long_name,
                         variable_standard_name,
                         variable_units,
                         variable_comment=None,
                         variable_type=float32,
                         variable_fill_value=9.96921e+36):

        var = self._create_nc_var(
            variable_name,
            variable_long_name,
            variable_standard_name,
            variable_units,
            variable_comment,
            variable_type=variable_type,
            variable_fill_value=variable_fill_value
        )
        if self.coverage.rank == 0:
            logging.info(f"[DefaultWriter] Writing variable '{variable_long_name}'")
        fn = getattr(self.coverage, function_name)
        self._gathering_var(fn(), var)

    def _write_time_2d_scalar(self,
                              function_name,
                              variable_name,
                              variable_long_name,
                              variable_standard_name,
                              variable_units,
                              variable_comment=None,
                              variable_type=float32,
                              variable_fill_value=9.96921e+36):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self._create_nc_var(
                variable_name,
                variable_long_name,
                variable_standard_name,
                variable_units,
                variable_comment,
                is_time_var=True,
                variable_type=variable_type,
                variable_fill_value=variable_fill_value)

            fn = getattr(self.coverage, function_name)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_long_name}' at time '{time}'")
                self._gathering_var(fn(time), var, time_index=time_index)
                time_index += 1
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def _write_time_2d_vector(self,
                              function_name,
                              variable_name,
                              eastward_variable_name,
                              eastward_variable_long_name,
                              eastward_variable_standard_name,
                              northward_variable_name,
                              northward_variable_long_name,
                              northward_variable_standard_name,
                              variable_units,
                              variable_comment=None,
                              variable_type=float32,
                              variable_fill_value=9.96921e+36
                              ):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            ucomp = self._create_nc_var(
                eastward_variable_name,
                eastward_variable_long_name,
                eastward_variable_standard_name,
                variable_units,
                variable_comment,
                is_time_var=True,
                variable_type=variable_type,
                variable_fill_value=variable_fill_value)

            vcomp = self._create_nc_var(
                northward_variable_name,
                northward_variable_long_name,
                northward_variable_standard_name,
                variable_units,
                variable_comment,
                is_time_var=True,
                variable_type=variable_type,
                variable_fill_value=variable_fill_value)

            fn = getattr(self.coverage, function_name)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_name}' at time '{time}'")

                local_data_u, local_data_v = fn(time)
                self._gathering_var(local_data_u, ucomp, time_index=time_index)
                self._gathering_var(local_data_v, vcomp, time_index=time_index)
                time_index += 1
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def _write_depth_2d_scalar(self,
                               function_name,
                               variable_name,
                               variable_long_name,
                               variable_standard_name,
                               variable_units,
                               variable_comment=None,
                               variable_type=float32,
                               variable_fill_value=9.96921e+36
                               ):

        if (isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self._create_nc_var(
                variable_name,
                variable_long_name,
                variable_standard_name,
                variable_units,
                variable_comment,
                is_depth_var=True,
                variable_type=variable_type,
                variable_fill_value=variable_fill_value)

            fn = getattr(self.coverage, function_name)

            level_index = 0
            for level in self.coverage.read_axis_z():
                self._gathering_var(fn(level), var, level_index=level_index)
                level_index += 1
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'LevelCoverage' or 'TimeLevelCoverage'")

    def _write_time_depth_2d_scalar(self,
                                    function_name,
                                    variable_name,
                                    variable_long_name,
                                    variable_standard_name,
                                    variable_units,
                                    variable_comment=None,
                                    variable_type=float32,
                                    variable_fill_value=9.96921e+36
                                    ):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            var = self._create_nc_var(
                variable_name,
                variable_long_name,
                variable_standard_name,
                variable_units,
                variable_comment,
                is_time_var=True,
                is_depth_var=True,
                variable_type=variable_type,
                variable_fill_value=variable_fill_value)

            fn = getattr(self.coverage, function_name)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_long_name}' at time '{time}'")
                level_index = 0
                for level in self.coverage.read_axis_z():
                    self._gathering_var(fn(time, level), var, time_index=time_index, level_index=level_index)
                    level_index += 1

                time_index += 1
        else:
            raise CoverageError("DefaultWriter", "The given coverage is not an instance of 'TimeLevelCoverage'")

    def _write_time_depth_2d_vector(self,
                                    function_name,
                                    variable_name,
                                    eastward_variable_name,
                                    eastward_variable_long_name,
                                    eastward_variable_standard_name,
                                    northward_variable_name,
                                    northward_variable_long_name,
                                    northward_variable_standard_name,
                                    variable_units,
                                    variable_comment=None,
                                    variable_type=float32,
                                    variable_fill_value=9.96921e+36
                                    ):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            ucomp = self._create_nc_var(
                eastward_variable_name,
                eastward_variable_long_name,
                eastward_variable_standard_name,
                variable_units,
                variable_comment,
                is_time_var=True,
                variable_type=variable_type,
                variable_fill_value=variable_fill_value)

            vcomp = self._create_nc_var(
                northward_variable_name,
                northward_variable_long_name,
                northward_variable_standard_name,
                variable_units,
                variable_comment,
                is_time_var=True,
                variable_type=variable_type,
                variable_fill_value=variable_fill_value)

            fn = getattr(self.coverage, function_name)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_name}' at time '{time}'")

                level_index = 0
                for level in self.coverage.read_axis_z():
                    local_data_u, local_data_v = fn(time, level)
                    self._gathering_var(local_data_u, ucomp, time_index=time_index, level_index=level_index)
                    self._gathering_var(local_data_v, vcomp, time_index=time_index, level_index=level_index)
                    level_index += 1
                time_index += 1
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeLevelCoverage'")

    # Variables
    def write_variable_mesh_size(self):

        self._write_2d_scalar(
            'read_variable_mesh_size',
            VariableDefinition.VARIABLE_NAME['mesh_size'],
            VariableDefinition.LONG_NAME['mesh_size'],
            VariableDefinition.STANDARD_NAME['mesh_size'],
            VariableDefinition.CANONICAL_UNITS['mesh_size']
        )

    def write_variable_mesh_size_factor(self):

        var = self._create_nc_var(VariableDefinition.VARIABLE_NAME['mesh_size_factor'],
                                  VariableDefinition.LONG_NAME['mesh_size_factor'],
                                  VariableDefinition.STANDARD_NAME['mesh_size_factor'],
                                  VariableDefinition.CANONICAL_UNITS['mesh_size_factor'])

        x = self.coverage.read_variable_x_mesh_size()
        y = self.coverage.read_variable_y_mesh_size()
        factor = np.divide(x, y)
        self._gathering_var(factor, var)

    def write_variable_2D_sea_binary_mask(self):

        self._write_2d_scalar(
            'read_variable_2D_sea_binary_mask',
            VariableDefinition.VARIABLE_NAME['2d_sea_binary_mask'],
            VariableDefinition.LONG_NAME['2d_sea_binary_mask'],
            VariableDefinition.STANDARD_NAME['2d_sea_binary_mask'],
            VariableDefinition.CANONICAL_UNITS['2d_sea_binary_mask'],
            variable_comment="1 = sea, 0 = land",
            variable_type=int16,
            variable_fill_value=-9999
        )

    def write_variable_wet_binary_mask(self):
        self._write_time_2d_scalar(
            'read_variable_2D_wet_binary_mask_at_time',
            VariableDefinition.VARIABLE_NAME['wet_binary_mask'],
            VariableDefinition.LONG_NAME['wet_binary_mask'],
            VariableDefinition.STANDARD_NAME['wet_binary_mask'],
            VariableDefinition.CANONICAL_UNITS['wet_binary_mask'],
            variable_comment="1 = sea, 0 = land",
            variable_type=int16,
            variable_fill_value=-9999
        )

    def write_variable_3D_sea_binary_mask(self):

        self._write_time_2d_scalar(
            'read_variable_2D_sea_binary_mask_at_time',
            VariableDefinition.VARIABLE_NAME['3d_sea_binary_mask'],
            VariableDefinition.LONG_NAME['3d_sea_binary_mask'],
            VariableDefinition.STANDARD_NAME['3d_sea_binary_mask'],
            VariableDefinition.CANONICAL_UNITS['3d_sea_binary_mask'],
            variable_comment="1 = sea, 0 = land",
            variable_type=int16,
            variable_fill_value=-9999
        )

    def write_variable_3D_land_binary_mask(self):

        self._write_time_2d_scalar(
            'read_variable_3D_land_binary_mask_at_time',
            VariableDefinition.VARIABLE_NAME['3d_land_binary_mask'],
            VariableDefinition.LONG_NAME['3d_land_binary_mask'],
            VariableDefinition.STANDARD_NAME['3d_land_binary_mask'],
            VariableDefinition.CANONICAL_UNITS['3d_land_binary_mask'],
            variable_comment="1 = sea, 0 = land",
            variable_type=int16,
            variable_fill_value=-9999
        )

    #################
    # HYDRO
    # 2D
    #################
    def write_variable_bathymetry(self):

        self._write_2d_scalar(
            'read_variable_bathymetry',
            VariableDefinition.VARIABLE_NAME['bathymetry'],
            VariableDefinition.LONG_NAME['bathymetry'],
            VariableDefinition.STANDARD_NAME['bathymetry'],
            VariableDefinition.CANONICAL_UNITS['bathymetry']
        )

    def write_variable_barotropic_sea_water_velocity(self):

        self._write_time_2d_vector(
            'read_variable_barotropic_sea_water_velocity_at_time',
            'Barotropic Sea Water Velocity',
            VariableDefinition.VARIABLE_NAME['barotropic_eastward_sea_water_velocity'],
            VariableDefinition.LONG_NAME['barotropic_eastward_sea_water_velocity'],
            VariableDefinition.STANDARD_NAME['barotropic_eastward_sea_water_velocity'],
            VariableDefinition.VARIABLE_NAME['barotropic_northward_sea_water_velocity'],
            VariableDefinition.LONG_NAME['barotropic_northward_sea_water_velocity'],
            VariableDefinition.STANDARD_NAME['barotropic_northward_sea_water_velocity'],
            VariableDefinition.CANONICAL_UNITS['barotropic_eastward_sea_water_velocity'],
            variable_comment="cur=sqrt(U**2+V**2)")

    def write_variable_barotropic_sea_water_speed(self):

        self._write_time_2d_scalar(
            'read_variable_barotropic_sea_water_speed_at_time',
            VariableDefinition.VARIABLE_NAME['barotropic_sea_water_speed'],
            VariableDefinition.LONG_NAME['barotropic_sea_water_speed'],
            VariableDefinition.STANDARD_NAME['barotropic_sea_water_speed'],
            VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_speed']
        )

    def write_variable_barotropic_sea_water_to_direction(self):

        self._write_time_2d_scalar(
            'read_variable_barotropic_sea_water_to_direction_at_time',
            VariableDefinition.VARIABLE_NAME['barotropic_sea_water_to_direction'],
            VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction'],
            VariableDefinition.STANDARD_NAME['barotropic_sea_water_to_direction'],
            VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_to_direction']
        )

    def write_variable_barotropic_sea_water_from_direction(self):

        self._write_time_2d_scalar(
            'read_variable_barotropic_sea_water_from_direction_at_time',
            VariableDefinition.VARIABLE_NAME['barotropic_sea_water_from_direction'],
            VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction'],
            VariableDefinition.STANDARD_NAME['barotropic_sea_water_from_direction'],
            VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_from_direction']
        )

    def write_variable_ocean_tracer_residence_time(self):
        self._write_time_2d_scalar(
            'read_variable_ocean_tracer_residence_time_at_time',
            VariableDefinition.VARIABLE_NAME['ocean_tracer_residence_time'],
            VariableDefinition.LONG_NAME['ocean_tracer_residence_time'],
            VariableDefinition.STANDARD_NAME['ocean_tracer_residence_time'],
            VariableDefinition.CANONICAL_UNITS['ocean_tracer_residence_time']
        )

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_surface_height_above_mean_sea_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_height_above_mean_sea_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_height_above_mean_sea_level'],
            VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'],
            VariableDefinition.STANDARD_NAME['sea_surface_height_above_mean_sea_level'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_mean_sea_level']
        )

    def write_variable_sea_surface_height_above_geoid(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_height_above_geoid_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_height_above_geoid'],
            VariableDefinition.LONG_NAME['sea_surface_height_above_geoid'],
            VariableDefinition.STANDARD_NAME['sea_surface_height_above_geoid'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_geoid']
        )

    def write_variable_sea_water_column_thickness(self):

        self._write_time_2d_scalar(
            'read_variable_sea_water_column_thickness_at_time',
            VariableDefinition.VARIABLE_NAME['sea_water_column_thickness'],
            VariableDefinition.LONG_NAME['sea_water_column_thickness'],
            VariableDefinition.STANDARD_NAME['sea_water_column_thickness'],
            VariableDefinition.CANONICAL_UNITS['sea_water_column_thickness']
        )

    def write_variable_sea_surface_temperature(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_temperature_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_temperature'],
            VariableDefinition.LONG_NAME['sea_surface_temperature'],
            VariableDefinition.STANDARD_NAME['sea_surface_temperature'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_temperature']
        )

    def write_variable_sea_surface_salinity(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_salinity_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_salinity'],
            VariableDefinition.LONG_NAME['sea_surface_salinity'],
            VariableDefinition.STANDARD_NAME['sea_surface_salinity'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_salinity']
        )

    def write_variable_sea_water_velocity_at_sea_water_surface(self):

        self._write_time_2d_vector(
            'read_variable_sea_water_velocity_at_sea_water_surface_at_time',
            'Sea Water Velocity at Sea Water Surface',
            VariableDefinition.VARIABLE_NAME['eastward_sea_water_velocity_at_sea_water_surface'],
            VariableDefinition.LONG_NAME['eastward_sea_water_velocity_at_sea_water_surface'],
            VariableDefinition.STANDARD_NAME['eastward_sea_water_velocity_at_sea_water_surface'],
            VariableDefinition.VARIABLE_NAME['northward_sea_water_velocity_at_sea_water_surface'],
            VariableDefinition.LONG_NAME['northward_sea_water_velocity_at_sea_water_surface'],
            VariableDefinition.STANDARD_NAME['northward_sea_water_velocity_at_sea_water_surface'],
            VariableDefinition.CANONICAL_UNITS['eastward_sea_water_velocity_at_sea_water_surface'],
            variable_comment="cur=sqrt(U**2+V**2)")

    def write_variable_sea_surface_residence_time(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_residence_time_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_residence_time'],
            VariableDefinition.LONG_NAME['sea_surface_residence_time'],
            VariableDefinition.STANDARD_NAME['sea_surface_residence_time'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_residence_time']
        )

    #################
    # HYDRO
    # Ground level
    #################

    def write_variable_sea_water_temperature_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_water_temperature_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_water_temperature_at_ground_level'],
            VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level'],
            VariableDefinition.STANDARD_NAME['sea_water_temperature_at_ground_level'],
            VariableDefinition.CANONICAL_UNITS['sea_water_temperature_at_ground_level']
        )

    def write_variable_sea_water_salinity_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_water_salinity_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_water_salinity_at_ground_level'],
            VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level'],
            VariableDefinition.STANDARD_NAME['sea_water_salinity_at_ground_level'],
            VariableDefinition.CANONICAL_UNITS['sea_water_salinity_at_ground_level']
        )

    def write_variable_sea_water_velocity_at_ground_level(self):

        self._write_time_2d_vector(
            'read_variable_sea_water_velocity_at_ground_level_at_time',
            'Sea Water Velocity at Ground Level',
            VariableDefinition.VARIABLE_NAME['eastward_sea_water_velocity_at_ground_level'],
            VariableDefinition.LONG_NAME['eastward_sea_water_velocity_at_ground_level'],
            VariableDefinition.STANDARD_NAME['eastward_sea_water_velocity_at_ground_level'],
            VariableDefinition.VARIABLE_NAME['northward_sea_water_velocity_at_ground_level'],
            VariableDefinition.LONG_NAME['northward_sea_water_velocity_at_ground_level'],
            VariableDefinition.STANDARD_NAME['northward_sea_water_velocity_at_ground_level'],
            VariableDefinition.CANONICAL_UNITS['eastward_sea_water_velocity_at_ground_level'],
            variable_comment="cur=sqrt(U**2+V**2)")

    def write_variable_residence_time_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_residence_time_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['residence_time_at_ground_level'],
            VariableDefinition.LONG_NAME['residence_time_at_ground_level'],
            VariableDefinition.STANDARD_NAME['residence_time_at_ground_level'],
            VariableDefinition.CANONICAL_UNITS['residence_time_at_ground_level']
        )

    #################
    # HYDRO
    # 3D
    #################

    def write_variable_depth(self):

        self._write_depth_2d_scalar(
            'read_variable_depth_at_depth',
            VariableDefinition.VARIABLE_NAME['depth_sigma'],
            VariableDefinition.LONG_NAME['depth_sigma'],
            VariableDefinition.STANDARD_NAME['depth_sigma'],
            VariableDefinition.CANONICAL_UNITS['depth_sigma']
        )

    def write_variable_sea_water_temperature(self):

        self._write_time_depth_2d_scalar(
            'read_variable_sea_water_temperature_at_time_and_depth',
            VariableDefinition.VARIABLE_NAME['sea_water_temperature'],
            VariableDefinition.LONG_NAME['sea_water_temperature'],
            VariableDefinition.STANDARD_NAME['sea_water_temperature'],
            VariableDefinition.CANONICAL_UNITS['sea_water_temperature']
        )

    def write_variable_sea_water_salinity(self):

        self._write_time_depth_2d_scalar(
            'read_variable_sea_water_salinity_at_time_and_depth',
            VariableDefinition.VARIABLE_NAME['sea_water_salinity'],
            VariableDefinition.LONG_NAME['sea_water_salinity'],
            VariableDefinition.STANDARD_NAME['sea_water_salinity'],
            VariableDefinition.CANONICAL_UNITS['sea_water_salinity']
        )

    def write_variable_baroclinic_sea_water_velocity(self):

        self._write_time_depth_2d_vector(
            'read_variable_baroclinic_sea_water_velocity_at_time_and_depth',
            'Baroclinic Sea Water Velocity',
            VariableDefinition.VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'],
            VariableDefinition.LONG_NAME['baroclinic_eastward_sea_water_velocity'],
            VariableDefinition.STANDARD_NAME['baroclinic_eastward_sea_water_velocity'],
            VariableDefinition.VARIABLE_NAME['baroclinic_northward_sea_water_velocity'],
            VariableDefinition.LONG_NAME['baroclinic_northward_sea_water_velocity'],
            VariableDefinition.STANDARD_NAME['baroclinic_northward_sea_water_velocity'],
            VariableDefinition.CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity'],
            variable_comment="cur=sqrt(U**2+V**2)")

    #################
    # WAVES
    # Sea Surface
    #################
    def write_variable_sea_surface_wave_significant_height(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_significant_height_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_significant_height'],
            VariableDefinition.LONG_NAME['sea_surface_wave_significant_height'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_significant_height'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_significant_height']
        )

    def write_variable_sea_surface_wave_breaking_height(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_breaking_height_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_breaking_height'],
            VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_breaking_height'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_breaking_height']
        )

    def write_variable_sea_surface_wave_mean_period(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_mean_period_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_mean_period'],
            VariableDefinition.LONG_NAME['sea_surface_wave_mean_period'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_mean_period'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_mean_period']
        )

    def write_variable_sea_surface_wave_peak_period(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_peak_period_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_peak_period'],
            VariableDefinition.LONG_NAME['sea_surface_wave_peak_period'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_peak_period'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_peak_period']
        )

    def write_variable_sea_surface_wave_from_direction(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_from_direction_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_from_direction'],
            VariableDefinition.LONG_NAME['sea_surface_wave_from_direction'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_from_direction'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_from_direction']
        )

    def write_variable_sea_surface_wave_to_direction(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_to_direction_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_to_direction'],
            VariableDefinition.LONG_NAME['sea_surface_wave_to_direction'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_to_direction'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_to_direction']
        )

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):

        self._write_time_2d_vector(
            'read_variable_sea_surface_wave_stokes_drift_velocity_at_time',
            'Surface Stokes Drift Velocity',
            VariableDefinition.VARIABLE_NAME['eastward_sea_surface_wave_stokes_drift_velocity'],
            VariableDefinition.LONG_NAME['eastward_sea_surface_wave_stokes_drift_velocity'],
            VariableDefinition.STANDARD_NAME['eastward_sea_surface_wave_stokes_drift_velocity'],
            VariableDefinition.VARIABLE_NAME['northward_sea_surface_wave_stokes_drift_velocity'],
            VariableDefinition.LONG_NAME['northward_sea_surface_wave_stokes_drift_velocity'],
            VariableDefinition.STANDARD_NAME['northward_sea_surface_wave_stokes_drift_velocity'],
            VariableDefinition.CANONICAL_UNITS['eastward_sea_surface_wave_stokes_drift_velocity'],
            variable_comment="cur=sqrt(U**2+V**2)")

    def write_variable_radiation_pressure_bernouilli_head(self):

        self._write_time_2d_scalar(
            'read_variable_radiation_pressure_bernouilli_head_at_time',
            VariableDefinition.VARIABLE_NAME['radiation_pressure_bernouilli_head'],
            VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head'],
            VariableDefinition.STANDARD_NAME['radiation_pressure_bernouilli_head'],
            VariableDefinition.CANONICAL_UNITS['radiation_pressure_bernouilli_head']
        )

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_energy_flux_to_ocean_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_flux_to_ocean'],
            VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_energy_flux_to_ocean'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_energy_flux_to_ocean']
        )

    #################
    # WAVES
    # Ground level
    #################
    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_dissipation_at_ground_level'],
            VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level'],
            VariableDefinition.STANDARD_NAME['sea_surface_wave_energy_dissipation_at_ground_level'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_wave_energy_dissipation_at_ground_level']
        )

    #################
    # WAVES
    # Momentum flux
    #################
    def write_variable_atmosphere_momentum_flux_to_waves(self):

        self._write_time_2d_vector(
            'read_variable_atmosphere_momentum_flux_to_waves_at_time',
            'Atmosphere Momentum Flux to Waves',
            VariableDefinition.VARIABLE_NAME['eastward_atmosphere_momentum_flux_to_waves'],
            VariableDefinition.LONG_NAME['eastward_atmosphere_momentum_flux_to_waves'],
            VariableDefinition.STANDARD_NAME['eastward_atmosphere_momentum_flux_to_waves'],
            VariableDefinition.VARIABLE_NAME['northward_atmosphere_momentum_flux_to_waves'],
            VariableDefinition.LONG_NAME['northward_atmosphere_momentum_flux_to_waves'],
            VariableDefinition.STANDARD_NAME['northward_atmosphere_momentum_flux_to_waves'],
            VariableDefinition.CANONICAL_UNITS['eastward_atmosphere_momentum_flux_to_waves'],
            variable_comment="cur=sqrt(U**2+V**2)")

    def write_variable_waves_momentum_flux_to_ocean(self):

        self._write_time_2d_vector(
            'read_variable_waves_momentum_flux_to_ocean_at_time',
            'Waves Momentum Flux To Ocean',
            VariableDefinition.VARIABLE_NAME['eastward_waves_momentum_flux_to_ocean'],
            VariableDefinition.LONG_NAME['eastward_waves_momentum_flux_to_ocean'],
            VariableDefinition.STANDARD_NAME['eastward_waves_momentum_flux_to_ocean'],
            VariableDefinition.VARIABLE_NAME['northward_waves_momentum_flux_to_ocean'],
            VariableDefinition.LONG_NAME['northward_waves_momentum_flux_to_ocean'],
            VariableDefinition.STANDARD_NAME['northward_waves_momentum_flux_to_ocean'],
            VariableDefinition.CANONICAL_UNITS['eastward_waves_momentum_flux_to_ocean'],
            variable_comment="cur=sqrt(U**2+V**2)")

    #################
    # METEO
    # 2D
    #################

    #################
    # METEO
    # Surface air
    #################
    def write_variable_sea_surface_air_pressure(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_air_pressure_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_air_pressure'],
            VariableDefinition.LONG_NAME['sea_surface_air_pressure'],
            VariableDefinition.STANDARD_NAME['sea_surface_air_pressure'],
            VariableDefinition.CANONICAL_UNITS['sea_surface_air_pressure']
        )

    def write_variable_surface_air_temperature(self):

        self._write_time_2d_scalar(
            'read_variable_surface_air_temperature_at_time',
            VariableDefinition.VARIABLE_NAME['surface_air_temperature'],
            VariableDefinition.LONG_NAME['surface_air_temperature'],
            VariableDefinition.STANDARD_NAME['surface_air_temperature'],
            VariableDefinition.CANONICAL_UNITS['surface_air_temperature']
        )

    def write_variable_dew_point_temperature(self):

        self._write_time_2d_scalar(
            'read_variable_dew_point_temperature_at_time',
            VariableDefinition.VARIABLE_NAME['dew_point_temperature'],
            VariableDefinition.LONG_NAME['dew_point_temperature'],
            VariableDefinition.STANDARD_NAME['dew_point_temperature'],
            VariableDefinition.CANONICAL_UNITS['dew_point_temperature']
        )

    def write_variable_rainfall_amount(self):

        self._write_time_2d_scalar(
            'read_variable_rainfall_amount_at_time',
            VariableDefinition.VARIABLE_NAME['rainfall_amount'],
            VariableDefinition.LONG_NAME['rainfall_amount'],
            VariableDefinition.STANDARD_NAME['rainfall_amount'],
            VariableDefinition.CANONICAL_UNITS['rainfall_amount']
        )

    def write_variable_surface_downward_sensible_heat_flux(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downward_sensible_heat_flux_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downward_sensible_heat_flux'],
            VariableDefinition.LONG_NAME['surface_downward_sensible_heat_flux'],
            VariableDefinition.STANDARD_NAME['surface_downward_sensible_heat_flux'],
            VariableDefinition.CANONICAL_UNITS['surface_downward_sensible_heat_flux']
        )

    def write_variable_surface_downward_latent_heat_flux(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downward_latent_heat_flux_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downward_latent_heat_flux'],
            VariableDefinition.LONG_NAME['surface_downward_latent_heat_flux'],
            VariableDefinition.STANDARD_NAME['surface_downward_latent_heat_flux'],
            VariableDefinition.CANONICAL_UNITS['surface_downward_latent_heat_flux']
        )

    def write_variable_surface_downward_solar_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downwards_solar_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downwards_solar_radiation'],
            VariableDefinition.LONG_NAME['surface_downwards_solar_radiation'],
            VariableDefinition.STANDARD_NAME['surface_downwards_solar_radiation'],
            VariableDefinition.CANONICAL_UNITS['surface_downwards_solar_radiation']
        )

    def write_variable_surface_downward_thermal_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downwards_thermal_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downwards_thermal_radiation'],
            VariableDefinition.LONG_NAME['surface_downwards_thermal_radiation'],
            VariableDefinition.STANDARD_NAME['surface_downwards_thermal_radiation'],
            VariableDefinition.CANONICAL_UNITS['surface_downwards_thermal_radiation']
        )

    def write_variable_surface_solar_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_solar_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_solar_radiation'],
            VariableDefinition.LONG_NAME['surface_solar_radiation'],
            VariableDefinition.STANDARD_NAME['surface_solar_radiation'],
            VariableDefinition.CANONICAL_UNITS['surface_solar_radiation']
        )

    def write_variable_surface_thermal_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_thermal_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_thermal_radiation'],
            VariableDefinition.LONG_NAME['surface_thermal_radiation'],
            VariableDefinition.STANDARD_NAME['surface_thermal_radiation'],
            VariableDefinition.CANONICAL_UNITS['surface_thermal_radiation']
        )

    def write_variable_wind_stress(self):

        self._write_time_2d_vector(
            'read_variable_wind_stress_at_time',
            'Wind Stress',
            VariableDefinition.VARIABLE_NAME['eastward_wind_stress'],
            VariableDefinition.LONG_NAME['eastward_wind_stress'],
            VariableDefinition.STANDARD_NAME['eastward_wind_stress'],
            VariableDefinition.VARIABLE_NAME['northward_wind_stress'],
            VariableDefinition.LONG_NAME['northward_wind_stress'],
            VariableDefinition.STANDARD_NAME['northward_wind_stress'],
            VariableDefinition.CANONICAL_UNITS['eastward_wind_stress'],
            variable_comment="cur=sqrt(U**2+V**2)")

    #################
    # METEO
    # At 10 m
    #################
    def write_variable_wind_10m(self):

        self._write_time_2d_vector(
            'read_variable_wind_10m_at_time',
            'Wind 10m',
            VariableDefinition.VARIABLE_NAME['eastward_wind_10m'],
            VariableDefinition.LONG_NAME['eastward_wind_10m'],
            VariableDefinition.STANDARD_NAME['eastward_wind_10m'],
            VariableDefinition.VARIABLE_NAME['northward_wind_10m'],
            VariableDefinition.LONG_NAME['northward_wind_10m'],
            VariableDefinition.STANDARD_NAME['northward_wind_10m'],
            VariableDefinition.CANONICAL_UNITS['eastward_wind_10m'],
            variable_comment="cur=sqrt(U**2+V**2)")

    def write_variable_wind_speed_10m(self):

        self._write_time_2d_scalar(
            'read_variable_wind_speed_10m_at_time',
            VariableDefinition.VARIABLE_NAME['wind_speed_10m'],
            VariableDefinition.LONG_NAME['wind_speed_10m'],
            VariableDefinition.STANDARD_NAME['wind_speed_10m'],
            VariableDefinition.CANONICAL_UNITS['wind_speed_10m']
        )

    def write_variable_wind_to_direction_10m(self):

        self._write_time_2d_scalar(
            'read_variable_wind_to_direction_10m_at_time',
            VariableDefinition.VARIABLE_NAME['wind_to_direction_10m'],
            VariableDefinition.LONG_NAME['wind_to_direction_10m'],
            VariableDefinition.STANDARD_NAME['wind_to_direction_10m'],
            VariableDefinition.CANONICAL_UNITS['wind_to_direction_10m']
        )

    def write_variable_wind_from_direction_10m(self):

        self._write_time_2d_scalar(
            'read_variable_wind_from_direction_10m_at_time',
            VariableDefinition.VARIABLE_NAME['wind_from_direction_10m'],
            VariableDefinition.LONG_NAME['wind_from_direction_10m'],
            VariableDefinition.STANDARD_NAME['wind_from_direction_10m'],
            VariableDefinition.CANONICAL_UNITS['wind_from_direction_10m']
        )
