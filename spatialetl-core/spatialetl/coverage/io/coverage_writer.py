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
from abc import ABC

import numpy as np


class CoverageWriter(ABC):

    def __init__(self, cov, myFile):
        self.coverage = cov;
        self.filename = myFile;

    def _gathering_var(self, local_data, var, time_index=None, level_index=None):
        if self.coverage.comm and self.coverage.rank != 0:
            self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # For MPI rank n°1
            if time_index is not None and level_index is not None:
                var[
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_t"].start + time_index:
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_y"],
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_x"]
                ] = local_data
            elif time_index is not None:
                var[
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_t"].start + time_index:
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_y"],
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_x"]
                ] = local_data
            elif level_index is not None:
                var[
                    level_index:level_index + 1,
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_y"],
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_x"]
                ] = local_data
            else:
                var[
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_y"],
                    self.coverage.parallel_map[self.coverage.rank]["dst_global_x"]
                ] = local_data

            if self.coverage.comm:
                # For others MPI rank
                for source in range(1, self.coverage.size):
                    recvbuf = np.empty([self.coverage.parallel_map[source]["dst_local_y_size"],
                                        self.coverage.parallel_map[source]["dst_local_x_size"]])
                    self.coverage.comm.Recv(recvbuf, source=source)

                    if time_index and level_index:
                        var[
                            self.coverage.parallel_map[source]["dst_global_t"].start + time_index:
                            self.coverage.parallel_map[source]["dst_global_t"].start + time_index + 1,
                            level_index:level_index + 1,
                            self.coverage.parallel_map[source]["dst_global_y"],
                            self.coverage.parallel_map[source]["dst_global_x"]
                        ] = recvbuf
                    elif time_index:
                        var[
                            self.coverage.parallel_map[source]["dst_global_t"].start + time_index:
                            self.coverage.parallel_map[source]["dst_global_t"].start + time_index + 1,
                            self.coverage.parallel_map[source]["dst_global_y"],
                            self.coverage.parallel_map[source]["dst_global_x"]
                        ] = recvbuf
                    elif level_index:
                        var[
                            level_index:level_index + 1,
                            self.coverage.parallel_map[source]["dst_global_y"],
                            self.coverage.parallel_map[source]["dst_global_x"]
                        ] = recvbuf
                    else:
                        var[
                            self.coverage.parallel_map[source]["dst_global_y"],
                            self.coverage.parallel_map[source]["dst_global_x"]
                        ] = recvbuf

    def close(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'close()'.")

    def write_variable_longitude(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_longitude()'.")

    def write_variable_latitude(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_latitude()'.")

    def write_variable_depth(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_depth()'.")

    def write_variable_time(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_time()'.")

    def write_variable_2D_sea_binary_mask(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_2D_sea_binary_mask()'.")

    def write_variable_2D_land_binary_mask(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_2D_land_binary_mask()'.")

    def write_variable_3D_sea_binary_mask(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_3D_sea_binary_mask()'.")

    def write_variable_3D_land_binary_mask(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_3D_land_binary_mask()'.")

    def write_variable_3D_land_binary_mask(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_3D_land_binary_mask()'.")

    def write_variable_wet_binary_mask(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wet_binary_mask()'.")

    def write_variable_mesh_size(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_mesh_size()'.")

    #################
    # HYDRO
    # 2D
    #################

    def write_variable_bathymetry(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_bathymetry()'.")

    def write_variable_barotropic_sea_water_velocity(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_barotropic_sea_water_velocity()'.")

    def write_variable_water_volume_transport_into_sea_water_from_rivers(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_water_volume_transport_into_sea_water_from_rivers()'.")

    #################
    # HYDRO
    # Sea Surface
    #################

    def write_variable_sea_surface_height_above_mean_sea_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_height_above_mean_sea_level()'.")

    def write_variable_sea_surface_height_above_geoid(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_height_above_geoid()'.")

    def write_variable_sea_surface_temperature(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_sea_surface_temperature()'.")

    def write_variable_sea_surface_salinity(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_sea_surface_salinity()'.")

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_pressure_at_sea_water_surface()'.")

    def write_variable_sea_surface_density(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_sea_surface_density()'.")

    def write_variable_sea_water_velocity_at_sea_water_surface(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_velocity_at_sea_water_surface()'.")

    #################
    # HYDRO
    # Ground level
    #################

    def write_variable_sea_water_temperature_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_temperature_at_ground_level()'.")

    def write_variable_sea_water_salinity_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_salinity_at_ground_level()'.")

    def write_variable_sea_water_pressure_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_pressure_at_ground_level()'.")

    def write_variable_sea_water_density_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_density_at_ground_level()'.")

    def write_variable_sea_water_velocity_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_velocity_at_ground_level()'.")

    #################
    # HYDRO
    # 3D
    #################

    def write_variable_sea_water_turbidity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_turbidity_and_depth()'.")

    def write_variable_sea_water_electrical_conductivity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_electrical_conductivity_and_depth()'.")

    def write_variable_sea_water_temperature_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_temperature_and_depth()'.")

    def write_variable_sea_water_salinity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_salinity_and_depth()'.")

    def write_variable_sea_water_density_and_depth(self, index_z):
        raise NotImplementedError(str(
            type(self)) + " don't have implemented the function 'write_variable_sea_water_density_and_depth()'.")

    def write_variable_baroclinic_sea_water_velocity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_baroclinic_sea_water_velocity_and_depth()'.")

    #################
    # WAVES
    # Sea Surface
    #################

    def write_variable_sea_surface_wave_significant_height(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_significant_height()'.")

    def write_variable_sea_surface_wave_breaking_height(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_breaking_height()'.")

    def write_variable_sea_surface_wave_mean_period(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_mean_period()'.")

    def write_variable_sea_surface_wave_peak_period(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_peak_period()'.")

    def write_variable_sea_surface_wave_from_direction(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_from_direction()'.")

    def write_variable_sea_surface_wave_to_direction(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_to_direction()'.")

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_stokes_drift_velocity()'.")

    def write_variable_radiation_pressure_bernouilli_head(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_radiation_pressure_bernouilli_head()'.")

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_energy_flux_to_ocean()'.")

    #################
    # WAVES
    # Ground level
    #################

    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_energy_dissipation_at_ground_level()'.")

    #################
    # WAVES
    # Momentum flux
    #################

    def write_variable_atmosphere_momentum_flux_to_waves(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_atmosphere_momentum_flux_to_waves()'.")

    def write_variable_waves_momentum_flux_to_ocean(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_waves_momentum_flux_to_ocean()'.")

    #################
    # METEO
    # 2D
    #################

    def write_variable_topography(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_topography()'.")

    def write_variable_rainfall_amount(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_rainfall_amount()'.")

    #################
    # METEO
    # Surface air
    #################

    def write_variable_surface_air_pressure(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_surface_air_pressure()'.")

    def write_variable_sea_surface_air_pressure(self):
        raise NotImplementedError(str(
            type(self)) + " don't have implemented the function 'write_variable_sea_surface_air_pressure()'.")

    def write_variable_wind_stress(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_stress()'.")

    def write_variable_surface_downward_sensible_heat_flux(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downward_sensible_heat_flux()'.")

    def write_variable_surface_downward_latent_heat_flux(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downward_latent_heat_flux()'.")

    def write_variable_surface_air_temperature(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_surface_air_temperature()'.")

    def write_variable_dew_point_temperature(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_dew_point_temperature()'.")

    def write_variable_surface_downward_solar_radiation(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downward_solar_radiation()'.")

    def write_variable_surface_downward_thermal_radiation(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downwards_thermal_radiation()'.")

    def write_variable_surface_solar_radiation(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_surface_solar_radiation()'.")

    def write_variable_surface_thermal_radiation(self):
        raise NotImplementedError(str(
            type(self)) + " don't have implemented the function 'write_variable_surface_thermal_radiation()'.")

    #################
    # METEO
    # At 10 m
    #################

    def write_variable_wind_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_10m()'.")

    def write_variable_wind_speed_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_speed_10m()'.")

    def write_variable_wind_to_direction_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_to_direction_10m()'.")

    def write_variable_wind_from_direction_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_from_direction_10m()'.")
