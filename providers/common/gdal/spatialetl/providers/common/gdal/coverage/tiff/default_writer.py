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

import os

import numpy as np
from osgeo import gdal
from osgeo import osr

from spatialetl.coverage import LevelCoverage
from spatialetl.coverage.io.coverage_writer import CoverageWriter
from spatialetl.coverage.time_coverage import TimeCoverage
from spatialetl.coverage.time_level_coverage import TimeLevelCoverage
from spatialetl.exception.coverage_error import CoverageError
from spatialetl.utils.logger import logging
from spatialetl.utils.variable_definition import VariableDefinition


class DefaultWriter(CoverageWriter):

    def __init__(self, cov, myFile):
        CoverageWriter.__init__(self, cov, myFile);

        if self.coverage.is_regular_grid() == False:
            raise ValueError("This writer supports only Coverage with a regular horizontal axis.")

        if os.path.isdir(self.filename) is False:
            raise ValueError("Filename has to be a directory.")

        gdal.AllRegister()
        self.driver = gdal.GetDriverByName('GTiff')

        self.rows = self.coverage.get_x_size(type="target_global")
        self.cols = self.coverage.get_y_size(type="target_global")

        xmin = np.min(self.coverage.read_axis_x(type="target_global"))
        xmax = np.max(self.coverage.read_axis_x(type="target_global"))
        ymin = np.min(self.coverage.read_axis_y(type="target_global"))
        ymax = np.max(self.coverage.read_axis_y(type="target_global"))
        x_pixel_size = round((xmax - xmin) / self.coverage.get_x_size(type="target_global"), 6)
        y_pixel_size = round((ymax - ymin) / self.coverage.get_y_size(type="target_global"), 6)

        self.geotransform = (xmin, x_pixel_size, 0, ymin, 0, y_pixel_size)

    def close(self):
        return

    def _create_file(self,
                     variable_name,
                     time=None,
                     depth=None,
                     variable_type=gdal.GDT_Float32,
                     variable_fill_value=9.96921e+36):
        if self.coverage.rank == 0:
            if time and depth:
                file = self.driver.Create(os.path.join(
                    self.filename,
                    f"{time.strftime('%Y%m%d_%H%M%S')}_{depth}m_{variable_name}.tiff"
                ),
                    int(self.rows), int(self.cols), 1, variable_type)
            elif depth:
                file = self.driver.Create(os.path.join(
                    self.filename,
                    f"{depth}m_{variable_name}.tiff"
                ),
                    int(self.rows), int(self.cols), 1, variable_type)
            elif time:
                file = self.driver.Create(os.path.join(
                    self.filename,
                    f"{time.strftime('%Y%m%d_%H%M%S')}_{variable_name}.tiff"
                ),
                    int(self.rows), int(self.cols), 1, variable_type)
            else:
                file = self.driver.Create(os.path.join(
                    self.filename,
                    f"{variable_name}.tiff"
                ),
                    int(self.rows), int(self.cols), 1, variable_type)
            # CRS info
            proj = osr.SpatialReference()
            proj.SetWellKnownGeogCS("EPSG:4326")
            file.SetProjection(proj.ExportToWkt())
            file.SetGeoTransform(self.geotransform)
            file.GetRasterBand(1).SetNoDataValue(variable_fill_value)
        else:
            file = None

        return file

    def _create_var(self,
                    is_time_var=False,
                    is_depth_var=False,
                    ):
        if self.coverage.rank == 0:
            if is_time_var and is_depth_var:
                var = np.empty(
                    [self.coverage.get_t_size(type="target_global"), self.coverage.get_z_size(type="target_global"),
                     int(self.cols), int(self.rows)])
                var[:] = np.nan
            elif is_depth_var:
                var = np.empty([self.coverage.get_z_size(type="target_global"), int(self.cols), int(self.rows)])
                var[:] = np.nan
            elif is_time_var:
                var = np.empty([self.coverage.get_t_size(type="target_global"), int(self.cols), int(self.rows)])
                var[:] = np.nan
            else:
                var = np.empty([int(self.cols), int(self.rows)])
                var[:] = np.nan
        else:
            var = None

        return var

    def _write_2d_scalar(self,
                         function_name,
                         variable_name,
                         variable_long_name,
                         variable_type=gdal.GDT_Float32,
                         variable_fill_value=9.96921e+36):

        file = self._create_file(
            variable_name,
            variable_type=variable_type,
            variable_fill_value=variable_fill_value
        )

        var = self._create_var()

        if self.coverage.rank == 0:
            logging.info(f"[DefaultWriter] Writing variable '{variable_long_name}'")
        fn = getattr(self.coverage, function_name)

        self._gathering_var(fn(), var)
        file.GetRasterBand(1).WriteArray(var)

    def _write_time_2d_scalar(self,
                              function_name,
                              variable_name,
                              variable_long_name,
                              variable_type=gdal.GDT_Float32,
                              variable_fill_value=9.96921e+36):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            fn = getattr(self.coverage, function_name)

            var = self._create_var(is_time_var=True)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_long_name}' at time '{time}'")

                file = self._create_file(
                    variable_name,
                    time=time,
                    variable_type=variable_type,
                    variable_fill_value=variable_fill_value)

                self._gathering_var(fn(time), var, time_index=time_index)
                file.GetRasterBand(1).WriteArray(var[time_index])
                time_index += 1
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def _write_time_2d_vector(self,
                              function_name,
                              variable_name,
                              eastward_variable_name,
                              northward_variable_name,
                              variable_type=gdal.GDT_Float32,
                              variable_fill_value=9.96921e+36
                              ):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            fn = getattr(self.coverage, function_name)

            ucomp = self._create_var(is_time_var=True)
            vcomp = self._create_var(is_time_var=True)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_name}' at time '{time}'")

                local_data_u, local_data_v = fn(time)

                ufile = self._create_file(
                    eastward_variable_name,
                    time=time,
                    variable_type=variable_type,
                    variable_fill_value=variable_fill_value)

                vfile = self._create_file(
                    northward_variable_name,
                    time=time,
                    variable_type=variable_type,
                    variable_fill_value=variable_fill_value)

                self._gathering_var(local_data_u, ucomp, time_index=time_index)
                self._gathering_var(local_data_v, vcomp, time_index=time_index)
                ufile.GetRasterBand(1).WriteArray(ucomp[time_index])
                vfile.GetRasterBand(1).WriteArray(vcomp[time_index])

                time_index += 1
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def _write_depth_2d_scalar(self,
                               function_name,
                               variable_name,
                               variable_long_name,
                               variable_type=gdal.GDT_Float32,
                               variable_fill_value=9.96921e+36
                               ):

        if (isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info(
                    f"[DefaultWriter] Writing variable '{variable_long_name}'")

            fn = getattr(self.coverage, function_name)

            var = self._create_var(is_depth_var=True)

            level_index = 0
            for level in self.coverage.read_axis_z():
                file = self._create_file(
                    variable_name,
                    depth=level,
                    variable_type=variable_type,
                    variable_fill_value=variable_fill_value)

                self._gathering_var(fn(level), var, level_index=level_index)
                file.GetRasterBand(1).WriteArray(var[level_index])

                level_index += 1
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'LevelCoverage' or 'TimeLevelCoverage'")

    def _write_time_depth_2d_scalar(self,
                                    function_name,
                                    variable_name,
                                    variable_long_name,
                                    variable_type=gdal.GDT_Float32,
                                    variable_fill_value=9.96921e+36
                                    ):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            fn = getattr(self.coverage, function_name)

            var = self._create_var(is_time_var=True, is_depth_var=True)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_long_name}' at time '{time}'")
                level_index = 0
                for level in self.coverage.read_axis_z():
                    file = self._create_file(
                        variable_name,
                        variable_long_name,
                        time=time,
                        depth=level,
                        variable_type=variable_type,
                        variable_fill_value=variable_fill_value)

                    self._gathering_var(fn(time, level), var, time_index=time_index, level_index=level_index)
                    file.GetRasterBand(1).WriteArray(var[time_index, level_index])

                    level_index += 1

                time_index += 1
        else:
            raise CoverageError("DefaultWriter", "The given coverage is not an instance of 'TimeLevelCoverage'")

    def _write_time_depth_2d_vector(self,
                                    function_name,
                                    variable_name,
                                    eastward_variable_name,
                                    northward_variable_name,
                                    variable_type=gdal.GDT_Float32,
                                    variable_fill_value=9.96921e+36
                                    ):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            fn = getattr(self.coverage, function_name)

            ucomp = self._create_var(is_time_var=True, is_depth_var=True)
            vcomp = self._create_var(is_time_var=True, is_depth_var=True)

            time_index = 0
            for time in self.coverage.read_axis_t():
                if self.coverage.rank == 0:
                    logging.info(
                        f"[DefaultWriter] Writing variable '{variable_name}' at time '{time}'")

                level_index = 0
                for level in self.coverage.read_axis_z():
                    local_data_u, local_data_v = fn(time, level)

                    ufile = self._create_file(
                        eastward_variable_name,
                        time=time,
                        variable_type=variable_type,
                        variable_fill_value=variable_fill_value)

                    vfile = self._create_file(
                        northward_variable_name,
                        time=time,
                        variable_type=variable_type,
                        variable_fill_value=variable_fill_value)

                    self._gathering_var(local_data_u, ucomp, time_index=time_index, level_index=level_index)
                    self._gathering_var(local_data_v, vcomp, time_index=time_index, level_index=level_index)
                    ufile.GetRasterBand(1).WriteArray(ucomp[time_index, level_index])
                    vfile.GetRasterBand(1).WriteArray(vcomp[time_index, level_index])

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
            VariableDefinition.LONG_NAME['mesh_size']
        )

    def write_variable_2D_sea_binary_mask(self):

        self._write_2d_scalar(
            'read_variable_2D_sea_binary_mask',
            VariableDefinition.VARIABLE_NAME['2d_sea_binary_mask'],
            VariableDefinition.LONG_NAME['2d_sea_binary_mask'],
            variable_type=gdal.GDT_Int16,
            variable_fill_value=-9999
        )

    def write_variable_wet_binary_mask(self):

        self._write_time_2d_scalar(
            'read_variable_2D_wet_binary_mask_at_time',
            VariableDefinition.VARIABLE_NAME['2d_sea_binary_mask'],
            VariableDefinition.LONG_NAME['2d_sea_binary_mask'],
            variable_type=gdal.GDT_Int16,
            variable_fill_value=-9999
        )

    def write_variable_3D_sea_binary_mask(self):

        self._write_time_2d_scalar(
            'read_variable_2D_sea_binary_mask_at_time',
            VariableDefinition.VARIABLE_NAME['3d_sea_binary_mask'],
            VariableDefinition.LONG_NAME['3d_sea_binary_mask'],
            variable_type=gdal.GDT_Int16,
            variable_fill_value=-9999
        )

    def write_variable_3D_land_binary_mask(self):

        self._write_time_2d_scalar(
            'read_variable_3D_land_binary_mask_at_time',
            VariableDefinition.VARIABLE_NAME['3d_land_binary_mask'],
            VariableDefinition.LONG_NAME['3d_land_binary_mask'],
            variable_type=gdal.GDT_Int16,
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
            VariableDefinition.LONG_NAME['bathymetry']
        )

    def write_variable_barotropic_sea_water_velocity(self):

        self._write_time_2d_vector(
            'read_variable_barotropic_sea_water_velocity_at_time',
            'Barotropic Sea Water Velocity',
            VariableDefinition.VARIABLE_NAME['barotropic_eastward_sea_water_velocity'],
            VariableDefinition.VARIABLE_NAME['barotropic_northward_sea_water_velocity']
        )

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
            VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction']
        )

    def write_variable_barotropic_sea_water_from_direction(self):
        self._write_time_2d_scalar(
            'read_variable_barotropic_sea_water_from_direction_at_time',
            VariableDefinition.VARIABLE_NAME['barotropic_sea_water_from_direction'],
            VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction']
        )

    def write_variable_ocean_tracer_residence_time(self):

        self._write_time_2d_scalar(
            'read_variable_ocean_tracer_residence_time_at_time',
            VariableDefinition.VARIABLE_NAME['ocean_tracer_residence_time'],
            VariableDefinition.LONG_NAME['ocean_tracer_residence_time']
        )

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_surface_height_above_mean_sea_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_height_above_mean_sea_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_height_above_mean_sea_level'],
            VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']
        )

    def write_variable_sea_surface_height_above_geoid(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_height_above_geoid_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_height_above_geoid'],
            VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']
        )

    def write_variable_sea_water_column_thickness(self):

        self._write_time_2d_scalar(
            'read_variable_sea_water_column_thickness_at_time',
            VariableDefinition.VARIABLE_NAME['sea_water_column_thickness'],
            VariableDefinition.LONG_NAME['sea_water_column_thickness']
        )

    def write_variable_sea_surface_temperature(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_temperature_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_temperature'],
            VariableDefinition.LONG_NAME['sea_surface_temperature']
        )

    def write_variable_sea_surface_salinity(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_salinity_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_salinity'],
            VariableDefinition.LONG_NAME['sea_surface_salinity']
        )

    def write_variable_sea_water_velocity_at_sea_water_surface(self):

        self._write_time_2d_vector(
            'read_variable_sea_water_velocity_at_sea_water_surface_at_time',
            'Sea Water Velocity at Sea Water Surface',
            VariableDefinition.VARIABLE_NAME['eastward_sea_water_velocity_at_sea_water_surface'],
            VariableDefinition.VARIABLE_NAME['northward_sea_water_velocity_at_sea_water_surface']
        )

    def write_variable_sea_surface_residence_time(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_residence_time_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_residence_time'],
            VariableDefinition.LONG_NAME['sea_surface_residence_time']
        )

    #################
    # HYDRO
    # Ground level
    #################

    def write_variable_sea_water_temperature_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_water_temperature_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_water_temperature_at_ground_level'],
            VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']
        )

    def write_variable_sea_water_salinity_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_water_salinity_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_water_salinity_at_ground_level'],
            VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']
        )

    def write_variable_sea_water_velocity_at_ground_level(self):

        self._write_time_2d_vector(
            'read_variable_sea_water_velocity_at_ground_level_at_time',
            'Sea Water Velocity at Ground Level',
            VariableDefinition.VARIABLE_NAME['eastward_sea_water_velocity_at_ground_level'],
            VariableDefinition.VARIABLE_NAME['northward_sea_water_velocity_at_ground_level']
        )

    def write_variable_residence_time_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_residence_time_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['residence_time_at_ground_level'],
            VariableDefinition.LONG_NAME['residence_time_at_ground_level']
        )

    #################
    # HYDRO
    # 3D
    #################

    def write_variable_depth(self):

        self._write_depth_2d_scalar(
            'read_variable_depth_at_depth',
            VariableDefinition.VARIABLE_NAME['depth_sigma'],
            VariableDefinition.LONG_NAME['depth_sigma']
        )

    def write_variable_sea_water_temperature(self):

        self._write_time_depth_2d_scalar(
            'read_variable_sea_water_temperature_at_time_and_depth',
            VariableDefinition.VARIABLE_NAME['sea_water_temperature'],
            VariableDefinition.LONG_NAME['sea_water_temperature'],
        )

    def write_variable_sea_water_salinity(self):

        self._write_time_depth_2d_scalar(
            'read_variable_sea_water_salinity_at_time_and_depth',
            VariableDefinition.VARIABLE_NAME['sea_water_salinity'],
            VariableDefinition.LONG_NAME['sea_water_salinity']
        )

    def write_variable_baroclinic_sea_water_velocity(self):

        self._write_time_depth_2d_vector(
            'read_variable_baroclinic_sea_water_velocity_at_time_and_depth',
            'Baroclinic Sea Water Velocity',
            VariableDefinition.VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'],
            VariableDefinition.VARIABLE_NAME['baroclinic_northward_sea_water_velocity'],
        )

    #################
    # WAVES
    # Sea Surface
    #################
    def write_variable_sea_surface_wave_significant_height(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_significant_height_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_significant_height'],
            VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']
        )

    def write_variable_sea_surface_wave_breaking_height(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_breaking_height_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_breaking_height'],
            VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']
        )

    def write_variable_sea_surface_wave_mean_period(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_mean_period_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_mean_period'],
            VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']
        )

    def write_variable_sea_surface_wave_peak_period(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_peak_period_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_peak_period'],
            VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']
        )

    def write_variable_sea_surface_wave_from_direction(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_from_direction_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_from_direction'],
            VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']
        )

    def write_variable_sea_surface_wave_to_direction(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_to_direction_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_to_direction'],
            VariableDefinition.LONG_NAME['sea_surface_wave_to_direction'],
        )

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):

        self._write_time_2d_vector(
            'read_variable_sea_surface_wave_stokes_drift_velocity_at_time',
            'Surface Stokes Drift Velocity',
            VariableDefinition.VARIABLE_NAME['eastward_sea_surface_wave_stokes_drift_velocity'],
            VariableDefinition.VARIABLE_NAME['northward_sea_surface_wave_stokes_drift_velocity'],
        )

    def write_variable_radiation_pressure_bernouilli_head(self):

        self._write_time_2d_scalar(
            'read_variable_radiation_pressure_bernouilli_head_at_time',
            VariableDefinition.VARIABLE_NAME['radiation_pressure_bernouilli_head'],
            VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']
        )

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_energy_flux_to_ocean_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_flux_to_ocean'],
            VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean'],
        )

    #################
    # WAVES
    # Ground level
    #################
    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):

        self._write_time_2d_scalar(
            'read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time',
            VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_dissipation_at_ground_level'],
            VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level']
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
            VariableDefinition.VARIABLE_NAME['northward_atmosphere_momentum_flux_to_waves'],
        )

    def write_variable_waves_momentum_flux_to_ocean(self):

        self._write_time_2d_vector(
            'read_variable_waves_momentum_flux_to_ocean_at_time',
            'Waves Momentum Flux To Ocean',
            VariableDefinition.VARIABLE_NAME['eastward_waves_momentum_flux_to_ocean'],
            VariableDefinition.VARIABLE_NAME['northward_waves_momentum_flux_to_ocean'],
        )

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
            VariableDefinition.LONG_NAME['sea_surface_air_pressure']
        )

    def write_variable_surface_air_temperature(self):

        self._write_time_2d_scalar(
            'read_variable_surface_air_temperature_at_time',
            VariableDefinition.VARIABLE_NAME['surface_air_temperature'],
            VariableDefinition.LONG_NAME['surface_air_temperature']
        )

    def write_variable_dew_point_temperature(self):

        self._write_time_2d_scalar(
            'read_variable_dew_point_temperature_at_time',
            VariableDefinition.VARIABLE_NAME['dew_point_temperature'],
            VariableDefinition.LONG_NAME['dew_point_temperature']
        )

    def write_variable_rainfall_amount(self):

        self._write_time_2d_scalar(
            'read_variable_rainfall_amount_at_time',
            VariableDefinition.VARIABLE_NAME['rainfall_amount'],
            VariableDefinition.LONG_NAME['rainfall_amount']
        )

    def write_variable_surface_downward_sensible_heat_flux(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downward_sensible_heat_flux_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downward_sensible_heat_flux'],
            VariableDefinition.LONG_NAME['surface_downward_sensible_heat_flux']
        )

    def write_variable_surface_downward_latent_heat_flux(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downward_latent_heat_flux_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downward_latent_heat_flux'],
            VariableDefinition.LONG_NAME['surface_downward_latent_heat_flux']
        )

    def write_variable_surface_downward_solar_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downwards_solar_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downwards_solar_radiation'],
            VariableDefinition.LONG_NAME['surface_downwards_solar_radiation']
        )

    def write_variable_surface_downward_thermal_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_downwards_thermal_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_downwards_thermal_radiation'],
            VariableDefinition.LONG_NAME['surface_downwards_thermal_radiation']
        )

    def write_variable_surface_solar_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_solar_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_solar_radiation'],
            VariableDefinition.LONG_NAME['surface_solar_radiation']
        )

    def write_variable_surface_thermal_radiation(self):

        self._write_time_2d_scalar(
            'read_variable_surface_thermal_radiation_at_time',
            VariableDefinition.VARIABLE_NAME['surface_thermal_radiation'],
            VariableDefinition.LONG_NAME['surface_thermal_radiation']
        )

    def write_variable_wind_stress(self):

        self._write_time_2d_vector(
            'read_variable_wind_stress_at_time',
            'Wind Stress',
            VariableDefinition.VARIABLE_NAME['eastward_wind_stress'],
            VariableDefinition.VARIABLE_NAME['northward_wind_stress'],
        )

    #################
    # METEO
    # At 10 m
    #################
    def write_variable_wind_10m(self):

        self._write_time_2d_vector(
            'read_variable_wind_10m_at_time',
            'Wind 10m',
            VariableDefinition.VARIABLE_NAME['eastward_wind_10m'],
            VariableDefinition.VARIABLE_NAME['northward_wind_10m'],
        )

    def write_variable_wind_speed_10m(self):

        self._write_time_2d_scalar(
            'read_variable_wind_speed_10m_at_time',
            VariableDefinition.VARIABLE_NAME['wind_speed_10m'],
            VariableDefinition.LONG_NAME['wind_speed_10m']
        )

    def write_variable_wind_to_direction_10m(self):

        self._write_time_2d_scalar(
            'read_variable_wind_to_direction_10m_at_time',
            VariableDefinition.VARIABLE_NAME['wind_to_direction_10m'],
            VariableDefinition.LONG_NAME['wind_to_direction_10m']
        )

    def write_variable_wind_from_direction_10m(self):

        self._write_time_2d_scalar(
            'read_variable_wind_from_direction_10m_at_time',
            VariableDefinition.VARIABLE_NAME['wind_from_direction_10m'],
            VariableDefinition.LONG_NAME['wind_from_direction_10m']
        )
