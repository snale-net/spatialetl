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
class VariableDefinition:

    STANDARD_NAME = {};
    LONG_NAME = {};
    VARIABLE_NAME = {}
    CANONICAL_UNITS = {};

    STANDARD_NAME['longitude'] = "longitude"
    LONG_NAME['longitude'] = "Longitude"
    VARIABLE_NAME['longitude'] = "longitude"
    CANONICAL_UNITS['longitude'] = "degree_east"

    STANDARD_NAME['latitude'] = "latitude"
    LONG_NAME['latitude'] = "Latitude"
    VARIABLE_NAME['latitude'] = "latitude"
    CANONICAL_UNITS['latitude'] = "degree_north"

    STANDARD_NAME['depth'] = "depth"
    LONG_NAME['depth'] = "Positive depth"
    VARIABLE_NAME['depth'] = "depth"
    CANONICAL_UNITS['depth'] = "m"

    STANDARD_NAME['depth_sigma'] = "depth"
    LONG_NAME['depth_sigma'] = "Sigma depth"
    VARIABLE_NAME['depth_sigma'] = "depth_sigma"
    CANONICAL_UNITS['depth_sigma'] = "m"

    STANDARD_NAME['time'] = "time"
    LONG_NAME['time'] = "Time UTC"
    VARIABLE_NAME['time'] = "time"
    CANONICAL_UNITS['time'] = "year-month-day hour:minute:second UTC"

    STANDARD_NAME['point'] = "point"
    LONG_NAME['point'] = "Geographical Point"
    VARIABLE_NAME['point'] = "point"
    CANONICAL_UNITS['point'] = "identifier"

    STANDARD_NAME['sea_binary_mask'] = "sea_binary_mask"
    LONG_NAME['sea_binary_mask'] = "Sea Binary Mask"
    VARIABLE_NAME['sea_binary_mask'] = "sea_mask"
    CANONICAL_UNITS['sea_binary_mask'] = "1"

    STANDARD_NAME['land_binary_mask'] = "land_binary_mask"
    LONG_NAME['land_binary_mask'] = "Land Binary Mask"
    VARIABLE_NAME['land_binary_mask'] = "land_mask"
    CANONICAL_UNITS['land_binary_mask'] = "1"

    STANDARD_NAME['2d_sea_binary_mask'] = "sea_binary_mask"
    LONG_NAME['2d_sea_binary_mask'] = "Sea Binary Mask"
    VARIABLE_NAME['2d_sea_binary_mask'] = "sea_mask"
    CANONICAL_UNITS['2d_sea_binary_mask'] = "1"

    STANDARD_NAME['2d_land_binary_mask'] = "land_binary_mask"
    LONG_NAME['2d_land_binary_mask'] = "Land Binary Mask"
    VARIABLE_NAME['2d_land_binary_mask'] = "land_mask"
    CANONICAL_UNITS['2d_land_binary_mask'] = "1"

    STANDARD_NAME['wet_binary_mask'] = "wet_binary_mask"
    LONG_NAME['wet_binary_mask'] = "Wet Binary Mask"
    VARIABLE_NAME['wet_binary_mask'] = "wet_mask"
    CANONICAL_UNITS['wet_binary_mask'] = "1"

    STANDARD_NAME['mesh_size'] = "mesh_size"
    LONG_NAME['mesh_size'] = "Mesh Size"
    VARIABLE_NAME['mesh_size'] = "mesh_size"
    CANONICAL_UNITS['mesh_size'] = "m"

    STANDARD_NAME['x_mesh_size'] = "x_mesh_size"
    LONG_NAME['x_mesh_size'] = "X Axis Mesh Size"
    VARIABLE_NAME['x_mesh_size'] = "x_mesh_size"
    CANONICAL_UNITS['x_mesh_size'] = "m"

    STANDARD_NAME['y_mesh_size'] = "y_mesh_size"
    LONG_NAME['y_mesh_size'] = "Y Axis Mesh Size"
    VARIABLE_NAME['y_mesh_size'] = "y_mesh_size"
    CANONICAL_UNITS['y_mesh_size'] = "m"

    STANDARD_NAME['mesh_size_factor'] = "mesh_size_factor"
    LONG_NAME['mesh_size_factor'] = "Mesh Size Factor"
    VARIABLE_NAME['mesh_size_factor'] = "mesh_size_factor"
    CANONICAL_UNITS['mesh_size_factor'] = "x/y"

    #################
    # HYDRO
    # Sea Surface
    #################
    STANDARD_NAME['sea_surface_height_above_mean_sea_level'] = "sea_surface_height_above_mean_sea_level"
    LONG_NAME['sea_surface_height_above_mean_sea_level'] = "Sea Surface Height Above Mean Sea Level"
    VARIABLE_NAME['sea_surface_height_above_mean_sea_level'] = "ssh_msl"
    CANONICAL_UNITS['sea_surface_height_above_mean_sea_level'] = "m"

    STANDARD_NAME['sea_surface_height_above_geoid'] = "sea_surface_height_above_geoid"
    LONG_NAME['sea_surface_height_above_geoid'] = "Sea Surface Height Above Geoid"
    VARIABLE_NAME['sea_surface_height_above_geoid'] = "ssh_geoid"
    CANONICAL_UNITS['sea_surface_height_above_geoid'] = "m"

    STANDARD_NAME['sea_water_column_thickness'] = "sea_water_column_thickness"
    LONG_NAME['sea_water_column_thickness'] = "Sea Water Column Thickness"
    VARIABLE_NAME['sea_water_column_thickness'] = "water_thickness"
    CANONICAL_UNITS['sea_water_column_thickness'] = "m"

    STANDARD_NAME['sea_surface_temperature'] = "sea_surface_temperature"
    LONG_NAME['sea_surface_temperature'] = "Sea Surface Temperature"
    VARIABLE_NAME['sea_surface_temperature'] = "sea_surface_temperature"
    CANONICAL_UNITS['sea_surface_temperature'] = "C"

    STANDARD_NAME['sea_surface_salinity'] = "sea_surface_salinity"
    LONG_NAME['sea_surface_salinity'] = "Sea Surface Salinity"
    VARIABLE_NAME['sea_surface_salinity'] = "sea_surface_salinity"
    CANONICAL_UNITS['sea_surface_salinity'] = "psu"

    STANDARD_NAME['sea_water_pressure_at_sea_water_surface'] = "sea_water_pressure_at_sea_water_surface"
    LONG_NAME['sea_water_pressure_at_sea_water_surface'] = "Sea Water Pressure At Sea Water Surface"
    VARIABLE_NAME['sea_water_pressure_at_sea_water_surface'] = "sea_surface_water_pressure"
    CANONICAL_UNITS['sea_water_pressure_at_sea_water_surface'] = "dbar"

    STANDARD_NAME['sea_surface_density'] = "sea_surface_density"
    LONG_NAME['sea_surface_density'] = "Sea Surface Density"
    VARIABLE_NAME['sea_surface_density'] = "sea_surface_density"
    CANONICAL_UNITS['sea_surface_density'] = "kg m-3"

    STANDARD_NAME['eastward_sea_water_velocity_at_sea_water_surface'] = "eastward_sea_water_velocity_at_sea_water_surface"
    LONG_NAME['eastward_sea_water_velocity_at_sea_water_surface'] = "Eastward Sea Water Velocity At Sea Water Surface"
    VARIABLE_NAME['eastward_sea_water_velocity_at_sea_water_surface'] = "u_sea_surface_vel"
    CANONICAL_UNITS['eastward_sea_water_velocity_at_sea_water_surface'] = "m s-1"

    STANDARD_NAME['northward_sea_water_velocity_at_sea_water_surface'] = "northward_sea_water_velocity_at_sea_water_surface"
    LONG_NAME['northward_sea_water_velocity_at_sea_water_surface'] = "Northward Sea Water Velocity At Sea Water Surface"
    VARIABLE_NAME['northward_sea_water_velocity_at_sea_water_surface'] = "v_sea_surface_vel"
    CANONICAL_UNITS['northward_sea_water_velocity_at_sea_water_surface'] = "m s-1"

    STANDARD_NAME['sea_water_speed_at_sea_water_surface'] = "sea_water_speed_at_sea_water_surface"
    LONG_NAME['sea_water_speed_at_sea_water_surface'] = "Sea Water Speed At Sea Water Surface"
    VARIABLE_NAME['sea_water_speed_at_sea_water_surface'] = "sea_surface_speed"
    CANONICAL_UNITS['sea_water_speed_at_sea_water_surface'] = "m s-1"

    STANDARD_NAME['sea_water_from_direction_at_sea_water_surface'] = "sea_water_from_direction_at_sea_water_surface"
    LONG_NAME['sea_water_from_direction_at_sea_water_surface'] = "Sea Water From Direction At Sea Water Surface"
    VARIABLE_NAME['sea_water_from_direction_at_sea_water_surface'] = "sea_surface_from_dir"
    CANONICAL_UNITS['sea_water_from_direction_at_sea_water_surface'] = "degree"  # from North=0"

    STANDARD_NAME['sea_water_to_direction_at_sea_water_surface'] = "sea_water_to_direction_at_sea_water_surface"
    LONG_NAME['sea_water_to_direction_at_sea_water_surface'] = "Sea Water To Direction At Sea Water Surface"
    VARIABLE_NAME['sea_water_to_direction_at_sea_water_surface'] = "sea_surface_to_dir"
    CANONICAL_UNITS['sea_water_to_direction_at_sea_water_surface'] = "degree"  # from North=0"

    STANDARD_NAME['sea_surface_residence_time'] = "sea_surface_residence_time"
    LONG_NAME['sea_surface_residence_time'] = "Sea Surface Residence Time"
    VARIABLE_NAME['sea_surface_residence_time'] = "sea_surface_residence_time"
    CANONICAL_UNITS['sea_surface_residence_time'] = "day"

    #################
    # HYDRO
    # Ground level
    #################
    STANDARD_NAME['sea_water_temperature_at_ground_level'] = "sea_water_temperature_at_ground_level"
    LONG_NAME['sea_water_temperature_at_ground_level'] = "Sea Water Temperature At Ground Level"
    VARIABLE_NAME['sea_water_temperature_at_ground_level'] = "sea_bottom_temperature"
    CANONICAL_UNITS['sea_water_temperature_at_ground_level'] = "degree"

    STANDARD_NAME['sea_water_salinity_at_ground_level'] = "sea_water_salinity_at_ground_level"
    LONG_NAME['sea_water_salinity_at_ground_level'] = "Sea Water Salinity At Ground Level"
    VARIABLE_NAME['sea_water_salinity_at_ground_level'] = "sea_bottom_salinity"
    CANONICAL_UNITS['sea_water_salinity_at_ground_level'] = "psu"

    STANDARD_NAME['sea_water_pressure_at_ground_level'] = "sea_water_pressure_at_ground_level"
    LONG_NAME['sea_water_pressure_at_ground_level'] = "Sea Water Pressure At Sea Water Surface"
    VARIABLE_NAME['sea_water_pressure_at_ground_level'] = "sea_bottom_pressure"
    CANONICAL_UNITS['sea_water_pressure_at_ground_level'] = "dbar"

    STANDARD_NAME['sea_water_density_at_ground_level'] = "sea_water_density_at_ground_level"
    LONG_NAME['sea_water_density_at_ground_level'] = "Sea Water Density At Ground Level"
    VARIABLE_NAME['sea_water_density_at_ground_level'] = "sea_water_density_at_ground_level"
    CANONICAL_UNITS['sea_water_density_at_ground_level'] = "kg m-3"

    STANDARD_NAME['eastward_sea_water_velocity_at_ground_level'] = "eastward_sea_water_velocity_at_ground_level"
    LONG_NAME['eastward_sea_water_velocity_at_ground_level'] = "Eastward Sea Water Velocity At Ground Level"
    VARIABLE_NAME['eastward_sea_water_velocity_at_ground_level'] = "u_sea_bottom_vel"
    CANONICAL_UNITS['eastward_sea_water_velocity_at_ground_level'] = "m s-1"

    STANDARD_NAME['northward_sea_water_velocity_at_ground_level'] = "northward_sea_water_velocity_at_ground_level"
    LONG_NAME['northward_sea_water_velocity_at_ground_level'] = "Northward Sea Water Velocity At Ground Level"
    VARIABLE_NAME['northward_sea_water_velocity_at_ground_level'] = "v_sea_bottom_vel"
    CANONICAL_UNITS['northward_sea_water_velocity_at_ground_level'] = "m s-1"

    STANDARD_NAME['sea_water_speed_at_ground_level'] = "sea_water_speed_at_ground_level"
    LONG_NAME['sea_water_speed_at_ground_level'] = "Sea Water Speed At Ground Level"
    VARIABLE_NAME['sea_water_speed_at_ground_level'] = "sea_bottom_speed"
    CANONICAL_UNITS['sea_water_speed_at_ground_level'] = "m s-1"

    STANDARD_NAME['sea_water_from_direction_at_ground_level'] = "sea_water_from_direction_at_ground_level"
    LONG_NAME['sea_water_from_direction_at_ground_level'] = "Sea Water From Direction At Ground Level"
    VARIABLE_NAME['sea_water_from_direction_at_ground_level'] = "sea_bottom_from_dir"
    CANONICAL_UNITS['sea_water_from_direction_at_ground_level'] = "degree"

    STANDARD_NAME['sea_water_to_direction_at_ground_level'] = "sea_water_to_direction_at_ground_level"
    LONG_NAME['sea_water_to_direction_at_ground_level'] = "Sea Water To Direction At Ground Level"
    VARIABLE_NAME['sea_water_to_direction_at_ground_level'] = "sea_bottom_to_dir"
    CANONICAL_UNITS['sea_water_to_direction_at_ground_level'] = "degree"

    STANDARD_NAME['residence_time_at_ground_level'] = "residence_time_at_ground_level"
    LONG_NAME['residence_time_at_ground_level'] = "Residence Time At Ground Level"
    VARIABLE_NAME['residence_time_at_ground_level'] = "residence_time_at_ground_level"
    CANONICAL_UNITS['residence_time_at_ground_level'] = "day"

    #################
    # HYDRO
    # 2D
    #################
    STANDARD_NAME['bathymetry'] = "bathymetry"
    LONG_NAME['bathymetry'] = "Bathymetry"
    VARIABLE_NAME['bathymetry'] = "bathymetry"
    CANONICAL_UNITS['bathymetry'] = "m"

    STANDARD_NAME['barotropic_eastward_sea_water_velocity'] = "barotropic_eastward_sea_water_velocity"
    LONG_NAME['barotropic_eastward_sea_water_velocity'] = "Barotropic Eastward Sea Water Velocity"
    VARIABLE_NAME['barotropic_eastward_sea_water_velocity'] = "u_sea_water_bar_vel"
    CANONICAL_UNITS['barotropic_eastward_sea_water_velocity'] = "m s-1"

    STANDARD_NAME['barotropic_northward_sea_water_velocity'] = "barotropic_northward_sea_water_velocity"
    LONG_NAME['barotropic_northward_sea_water_velocity'] = "Barotropic Northward Sea Water Velocity"
    VARIABLE_NAME['barotropic_northward_sea_water_velocity'] = "v_sea_water_bar_vel"
    CANONICAL_UNITS['barotropic_northward_sea_water_velocity'] = "m s-1"

    STANDARD_NAME['barotropic_sea_water_speed'] = "barotropic_sea_water_speed"
    LONG_NAME['barotropic_sea_water_speed'] = "Barotropic Sea Water Speed"
    VARIABLE_NAME['barotropic_sea_water_speed'] = "sea_water_bar_speed"
    CANONICAL_UNITS['barotropic_sea_water_speed'] = "m s-1"

    STANDARD_NAME['barotropic_sea_water_from_direction'] = "barotropic_sea_water_from_direction"
    LONG_NAME['barotropic_sea_water_from_direction'] = "Barotropic Sea Water From Direction"
    VARIABLE_NAME['barotropic_sea_water_from_direction'] = "sea_water_bar_from_dir"
    CANONICAL_UNITS['barotropic_sea_water_from_direction'] = "degree"  # from North=0"

    STANDARD_NAME['barotropic_sea_water_to_direction'] = "barotropic_sea_water_to_direction"
    LONG_NAME['barotropic_sea_water_to_direction'] = "Barotropic Sea Water To Direction"
    VARIABLE_NAME['barotropic_sea_water_to_direction'] = "sea_water_bar_to_dir"
    CANONICAL_UNITS['barotropic_sea_water_to_direction'] = "degree"  # from North=0"

    STANDARD_NAME['water_volume_transport_into_sea_water_from_rivers'] = "water_volume_transport_into_sea_water_from_rivers"
    LONG_NAME['water_volume_transport_into_sea_water_from_rivers'] = "Water Volume Transport Into Sea Water From Rivers"
    VARIABLE_NAME['water_volume_transport_into_sea_water_from_rivers'] = "rivers_flux"
    CANONICAL_UNITS['water_volume_transport_into_sea_water_from_rivers'] = "m3 s-1"

    STANDARD_NAME['ocean_tracer_residence_time'] = "ocean_tracer_residence_time"
    LONG_NAME['ocean_tracer_residence_time'] = "Water Residence Time From Ocean Tracers"
    VARIABLE_NAME['ocean_tracer_residence_time'] = "residence_time"
    CANONICAL_UNITS['ocean_tracer_residence_time'] = "days"

    #################
    # HYDRO
    # 3D
    #################
    STANDARD_NAME['sea_water_turbidity'] = "sea_water_turbidity"
    LONG_NAME['sea_water_turbidity'] = "Sea Water Turbidity"
    VARIABLE_NAME['sea_water_turbidity'] = "sea_water_turbidity"
    CANONICAL_UNITS['sea_water_turbidity'] = "FTU"

    STANDARD_NAME['sea_water_electrical_conductivity'] = "sea_water_electrical_conductivity"
    LONG_NAME['sea_water_electrical_conductivity'] = "Sea Water Electrical Conductivity"
    VARIABLE_NAME['sea_water_electrical_conductivity'] = "sea_water_electrical_conductivity"
    CANONICAL_UNITS['sea_water_electrical_conductivity'] = "S m-1"

    STANDARD_NAME['sea_water_temperature'] = "sea_water_temperature"
    LONG_NAME['sea_water_temperature'] = "Sea Water Temperature"
    VARIABLE_NAME['sea_water_temperature'] = "sea_water_temperature"
    CANONICAL_UNITS['sea_water_temperature'] = "degree"

    STANDARD_NAME['sea_water_salinity'] = "sea_water_salinity"
    LONG_NAME['sea_water_salinity'] = "Sea Water Salinity"
    VARIABLE_NAME['sea_water_salinity'] = "sea_water_salinity"
    CANONICAL_UNITS['sea_water_salinity'] = "psu"

    STANDARD_NAME['sea_water_density'] = "sea_water_density"
    LONG_NAME['sea_water_density'] = "Sea Water Density"
    VARIABLE_NAME['sea_water_density'] = "sea_water_density"
    CANONICAL_UNITS['sea_water_density'] = "kg m-3"

    STANDARD_NAME['baroclinic_eastward_sea_water_velocity'] = "baroclinic_eastward_sea_water_velocity"
    LONG_NAME['baroclinic_eastward_sea_water_velocity'] = "Baroclinic Eastward Sea Water Velocity"
    VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'] = "u_sea_water_vel"
    CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity'] = "m s-1"

    STANDARD_NAME['baroclinic_northward_sea_water_velocity'] = "baroclinic_northward_sea_water_velocity"
    LONG_NAME['baroclinic_northward_sea_water_velocity'] = "Baroclinic Northward Sea Water Velocity"
    VARIABLE_NAME['baroclinic_northward_sea_water_velocity'] = "v_sea_water_vel"
    CANONICAL_UNITS['baroclinic_northward_sea_water_velocity'] = "m s-1"

    STANDARD_NAME['baroclinic_sea_water_speed'] = "baroclinic_sea_water_speed"
    LONG_NAME['baroclinic_sea_water_speed'] = "Baroclinic Sea Water Speed"
    VARIABLE_NAME['baroclinic_sea_water_speed'] = "sea_water_speed"
    CANONICAL_UNITS['baroclinic_sea_water_speed'] = "m s-1"

    STANDARD_NAME['baroclinic_sea_water_from_direction'] = "baroclinic_sea_water_from_direction"
    LONG_NAME['baroclinic_sea_water_from_direction'] = "Baroclinic Sea Water From Direction"
    VARIABLE_NAME['baroclinic_sea_water_from_direction'] = "sea_water_from_dir"
    CANONICAL_UNITS['baroclinic_sea_water_from_direction'] = "degree"  # from North=0"

    STANDARD_NAME['baroclinic_sea_water_to_direction'] = "baroclinic_sea_water_to_direction"
    LONG_NAME['baroclinic_sea_water_to_direction'] = "Baroclinic Sea Water To Direction"
    VARIABLE_NAME['baroclinic_sea_water_to_direction'] = "sea_water_to_dir"
    CANONICAL_UNITS['baroclinic_sea_water_to_direction'] = "degree"  # from North=0"

    #################
    # WAVES
    # Sea Surface
    #################
    STANDARD_NAME['sea_surface_wave_significant_height'] = "sea_surface_wave_significant_height"
    LONG_NAME['sea_surface_wave_significant_height'] = "Sea Surface Wave Significant Height"
    VARIABLE_NAME['sea_surface_wave_significant_height'] = "hs"
    CANONICAL_UNITS['sea_surface_wave_significant_height'] = "m"

    STANDARD_NAME['sea_surface_wave_breaking_height'] = "sea_surface_wave_breaking_height"
    LONG_NAME['sea_surface_wave_breaking_height'] = "Sea Surface Wave Breaking Height"
    VARIABLE_NAME['sea_surface_wave_breaking_height'] = "wch"
    CANONICAL_UNITS['sea_surface_wave_breaking_height'] = "m"

    STANDARD_NAME['sea_surface_wave_mean_period'] = "sea_surface_wave_mean_period"
    LONG_NAME['sea_surface_wave_mean_period'] = "Sea Surface Wave Mean Period"
    VARIABLE_NAME['sea_surface_wave_mean_period'] = "wave_mean_period"
    CANONICAL_UNITS['sea_surface_wave_mean_period'] = "s"

    STANDARD_NAME['sea_surface_wave_peak_period'] = "sea_surface_wave_peak_period"
    LONG_NAME['sea_surface_wave_peak_period'] = "Sea Surface Wave Peak Period"
    VARIABLE_NAME['sea_surface_wave_peak_period'] = "wave_peak_period"
    CANONICAL_UNITS['sea_surface_wave_peak_period'] = "s"

    STANDARD_NAME['sea_surface_wave_from_direction'] = "sea_surface_wave_from_direction"
    LONG_NAME['sea_surface_wave_from_direction'] = "Sea Surface Wave From Direction"
    VARIABLE_NAME['sea_surface_wave_from_direction'] = "wave_from_dir"
    CANONICAL_UNITS['sea_surface_wave_from_direction'] = "degree" # from North=0 / East=90"

    STANDARD_NAME['sea_surface_wave_to_direction'] = "sea_surface_wave_to_direction"
    LONG_NAME['sea_surface_wave_to_direction'] = "Sea Surface Wave To Direction"
    VARIABLE_NAME['sea_surface_wave_to_direction'] = "wave_to_dir"
    CANONICAL_UNITS['sea_surface_wave_to_direction'] = "degree" # from North=0 / East=90"

    STANDARD_NAME['eastward_sea_surface_wave_stokes_drift_velocity'] = "eastward_sea_surface_wave_stokes_drift_velocity"
    LONG_NAME['eastward_sea_surface_wave_stokes_drift_velocity'] = "Eastward Sea Surface Wave Stokes Drift Velocity"
    VARIABLE_NAME['eastward_sea_surface_wave_stokes_drift_velocity'] = "u_stokes_drift_vel"
    CANONICAL_UNITS['eastward_sea_surface_wave_stokes_drift_velocity'] = "m s-1"

    STANDARD_NAME['northward_sea_surface_wave_stokes_drift_velocity'] = "northward_sea_surface_wave_stokes_drift_velocity"
    LONG_NAME['northward_sea_surface_wave_stokes_drift_velocity'] = "Northward Sea Surface Wave Stokes Drift Velocity"
    VARIABLE_NAME['northward_sea_surface_wave_stokes_drift_velocity'] = "v_stokes_drift_vel"
    CANONICAL_UNITS['northward_sea_surface_wave_stokes_drift_velocity'] = "m s-1"

    STANDARD_NAME['radiation_pressure_bernouilli_head'] = "radiation_pressure_bernouilli_head"
    LONG_NAME['radiation_pressure_bernouilli_head'] = "Radiation Pressure Bernouilli_Head"
    VARIABLE_NAME['radiation_pressure_bernouilli_head'] = "pres_bernouilli_head"
    CANONICAL_UNITS['radiation_pressure_bernouilli_head'] = "m2 s-2"

    STANDARD_NAME['sea_surface_wave_energy_flux_to_ocean'] = "sea_surface_wave_energy_flux_to_ocean"
    LONG_NAME['sea_surface_wave_energy_flux_to_ocean'] = "Sea Surface Wave Energy Flux To Ocean"
    VARIABLE_NAME['sea_surface_wave_energy_flux_to_ocean'] = "wave_energy_flux_to_ocean"
    CANONICAL_UNITS['sea_surface_wave_energy_flux_to_ocean'] = "W m-2"

    STANDARD_NAME['sea_surface_wave_energy_dissipation_at_ground_level'] = "sea_surface_wave_energy_dissipation_at_ground_level"
    LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level'] = "Sea Surface Wave Energy Dissipation At Ground Level"
    VARIABLE_NAME['sea_surface_wave_energy_dissipation_at_ground_level'] = "wave_energy_dissip_bottom"
    CANONICAL_UNITS['sea_surface_wave_energy_dissipation_at_ground_level'] = "W m-2"

    #################
    # WAVES
    # Momentum flux
    #################
    STANDARD_NAME['eastward_atmosphere_momentum_flux_to_waves'] = "eastward_atmosphere_momentum_flux_to_waves"
    LONG_NAME['eastward_atmosphere_momentum_flux_to_waves'] = "Eastward Atmosphere Momentum Flux To Waves"
    VARIABLE_NAME['eastward_atmosphere_momentum_flux_to_waves'] = "u_taw"
    CANONICAL_UNITS['eastward_atmosphere_momentum_flux_to_waves'] = "m2 s-1"

    STANDARD_NAME['northward_atmosphere_momentum_flux_to_waves'] = "northward_atmosphere_momentum_flux_to_waves"
    LONG_NAME['northward_atmosphere_momentum_flux_to_waves'] = "Northward Atmosphere Momentum Flux To Waves"
    VARIABLE_NAME['northward_atmosphere_momentum_flux_to_waves'] = "v_taw"
    CANONICAL_UNITS['northward_atmosphere_momentum_flux_to_waves'] = "m2 s-1"

    STANDARD_NAME['eastward_waves_momentum_flux_to_ocean'] = "eastward_waves_momentum_flux_to_ocean"
    LONG_NAME['eastward_waves_momentum_flux_to_ocean'] = "Eastward Waves Momentum Flux To Ocean"
    VARIABLE_NAME['eastward_waves_momentum_flux_to_ocean'] = "u_two"
    CANONICAL_UNITS['eastward_waves_momentum_flux_to_ocean'] = "m2 s-1"

    STANDARD_NAME['northward_waves_momentum_flux_to_ocean'] = "northward_waves_momentum_flux_to_ocean"
    LONG_NAME['northward_waves_momentum_flux_to_ocean'] = "Northward Waves Momentum Flux To Ocean"
    VARIABLE_NAME['northward_waves_momentum_flux_to_ocean'] = "v_two"
    CANONICAL_UNITS['northward_waves_momentum_flux_to_ocean'] = "m2 s-1"

    #################
    # WAVES
    # Interaction with coastal structures
    #################
    STANDARD_NAME['wave_collision'] = "wave_collision"
    LONG_NAME['wave_collision'] = "Wave collision"
    VARIABLE_NAME['wave_collision'] = "wave_collision"
    CANONICAL_UNITS['wave_collision'] = "num h-1"

    STANDARD_NAME['wave_overtopping'] = "wave_overtopping"
    LONG_NAME['wave_overtopping'] = "Wave overtopping"
    VARIABLE_NAME['wave_overtopping'] = "wave_overtopping"
    CANONICAL_UNITS['wave_overtopping'] = "num h-1"

    #################
    # METEO
    # 2D
    #################
    STANDARD_NAME['topography'] = "topography"
    LONG_NAME['topography'] = "Topography"
    VARIABLE_NAME['topography'] = "topography"
    CANONICAL_UNITS['topography'] = "m"

    STANDARD_NAME['rainfall_amount'] = "rainfall_amount"
    LONG_NAME['rainfall_amount'] = "Rainfall Amount"
    VARIABLE_NAME['rainfall_amount'] = "rainfall"
    CANONICAL_UNITS['rainfall_amount'] = "kg m-1"

    #################
    # METEO
    # Surface air
    #################
    STANDARD_NAME['surface_air_pressure'] = "surface_air_pressure"
    LONG_NAME['surface_air_pressure'] = "Surface Air Pressure"
    VARIABLE_NAME['surface_air_pressure'] = "surface_air_pressure"
    CANONICAL_UNITS['surface_air_pressure'] = "hPa"

    STANDARD_NAME['sea_surface_air_pressure'] = "sea_surface_air_pressure"
    LONG_NAME['sea_surface_air_pressure'] = "Sea Surface Air Pressure"
    VARIABLE_NAME['sea_surface_air_pressure'] = "sea_surface_air_pressure"
    CANONICAL_UNITS['sea_surface_air_pressure'] = "hPa"

    STANDARD_NAME['eastward_wind_stress'] = "eastward_wind_stress"
    LONG_NAME['eastward_wind_stress'] = "Eastward Wind Stress"
    VARIABLE_NAME['eastward_wind_stress'] = "u_wind_stress"
    CANONICAL_UNITS['eastward_wind_stress'] = "N m2"

    STANDARD_NAME['northward_wind_stress'] = "northward_wind_stress"
    LONG_NAME['northward_wind_stress'] = "Northward Wind Stress"
    VARIABLE_NAME['northward_wind_stress'] = "v_wind_stress"
    CANONICAL_UNITS['northward_wind_stress'] = "N m2"

    STANDARD_NAME['wind_stress_stress'] = "wind_stress_stress"
    LONG_NAME['wind_stress_stress'] = "Wind Stress"
    VARIABLE_NAME['wind_stress_stress'] = "wind_stress_stress"
    CANONICAL_UNITS['wind_stress_stress'] = "N m2"

    STANDARD_NAME['wind_stress_to_direction'] = "wind_stress_to_direction"
    LONG_NAME['wind_stress_to_direction'] = "Wind Stress to Direction"
    VARIABLE_NAME['wind_stress_to_direction'] = "wind_stress_to_dir"
    CANONICAL_UNITS['wind_stress_to_direction'] = "degree"  # to North=0 / East=90"

    STANDARD_NAME['wind_stress_from_direction'] = "wind_stress_from_direction"
    LONG_NAME['wind_stress_from_direction'] = "Wind Stress From Direction"
    VARIABLE_NAME['wind_stress_from_direction'] = "wind_stress_from_dir"
    CANONICAL_UNITS['wind_stress_from_direction'] = "degree"  # from North=0 / East=90"

    STANDARD_NAME['surface_downward_sensible_heat_flux'] = "surface_downward_sensible_heat_flux"
    LONG_NAME['surface_downward_sensible_heat_flux'] = "Surface Downward Sensible Heat Flux"
    VARIABLE_NAME['surface_downward_sensible_heat_flux'] = "surface_downward_sensible_heat_flux"
    CANONICAL_UNITS['surface_downward_sensible_heat_flux'] = "W m2"

    STANDARD_NAME['surface_downward_latent_heat_flux'] = "surface_downward_latent_heat_flux"
    LONG_NAME['surface_downward_latent_heat_flux'] = "Surface Downward Latent Heat Flux"
    VARIABLE_NAME['surface_downward_latent_heat_flux'] = "surface_downward_latent_heat_flux"
    CANONICAL_UNITS['surface_downward_latent_heat_flux'] = "W m2"

    STANDARD_NAME['surface_air_temperature'] = "surface_air_temperature"
    LONG_NAME['surface_air_temperature'] = "Surface Air Temperature"
    VARIABLE_NAME['surface_air_temperature'] = "surface_air_temperature"
    CANONICAL_UNITS['surface_air_temperature'] = "C"

    STANDARD_NAME['dew_point_temperature'] = "dew_point_temperature"
    LONG_NAME['dew_point_temperature'] = "Dew Point Temperature"
    VARIABLE_NAME['dew_point_temperature'] = "dew_point_temperature"
    CANONICAL_UNITS['dew_point_temperature'] = "C"

    STANDARD_NAME['surface_downward_solar_radiation'] = "surface_downward_solar_radiation"
    LONG_NAME['surface_downward_solar_radiation'] = "Surface Downward Solar Radiation"
    VARIABLE_NAME['surface_downward_solar_radiation'] = "surface_downward_solar_radiation"
    CANONICAL_UNITS['surface_downward_solar_radiation'] = "W m2"

    STANDARD_NAME['surface_downward_thermal_radiation'] = "surface_downward_thermal_radiation"
    LONG_NAME['surface_downward_thermal_radiation'] = "Surface Downward Thermal Radiation"
    VARIABLE_NAME['surface_downward_thermal_radiation'] = "surface_downward_thermal_radiation"
    CANONICAL_UNITS['surface_downward_thermal_radiation'] = "W m2"

    STANDARD_NAME['surface_solar_radiation'] = "surface_solar_radiation"
    LONG_NAME['surface_solar_radiation'] = "Surface Solar Radiation"
    VARIABLE_NAME['surface_solar_radiation'] = "surface_solar_radiation"
    CANONICAL_UNITS['surface_solar_radiation'] = "W m2"

    STANDARD_NAME['surface_thermal_radiation'] = "surface_thermal_radiation"
    LONG_NAME['surface_thermal_radiation'] = "Surface Thermal Radiation"
    VARIABLE_NAME['surface_thermal_radiation'] = "surface_thermal_radiation"
    CANONICAL_UNITS['surface_thermal_radiation'] = "W m2"

    #################
    # METEO
    # At 10 m
    #################
    STANDARD_NAME['eastward_wind_10m'] = "eastward_wind_10m"
    LONG_NAME['eastward_wind_10m'] = "Eastward Wind 10m"
    VARIABLE_NAME['eastward_wind_10m'] = "u_wind_10m"
    CANONICAL_UNITS['eastward_wind_10m'] = "m s-1"

    STANDARD_NAME['northward_wind_10m'] = "northward_wind_10m"
    LONG_NAME['northward_wind_10m'] = "Northward Wind 10m"
    VARIABLE_NAME['northward_wind_10m'] = "v_wind_10m"
    CANONICAL_UNITS['northward_wind_10m'] = "m s-1"

    STANDARD_NAME['wind_speed_10m'] = "wind_speed_10m"
    LONG_NAME['wind_speed_10m'] = "Wind Speed 10m"
    VARIABLE_NAME['wind_speed_10m'] = "wind_speed_10m"
    CANONICAL_UNITS['wind_speed_10m'] = "m s-1"

    STANDARD_NAME['wind_to_direction_10m'] = "wind_to_direction_10m"
    LONG_NAME['wind_to_direction_10m'] = "Wind to Direction 10m"
    VARIABLE_NAME['wind_to_direction_10m'] = "wind_to_dir_10m"
    CANONICAL_UNITS['wind_to_direction_10m'] = "degree" # to North=0 / East=90"

    STANDARD_NAME['wind_from_direction_10m'] = "wind_from_direction_10m"
    LONG_NAME['wind_from_direction_10m'] = "Wind From Direction 10m"
    VARIABLE_NAME['wind_from_direction_10m'] = "wind_from_dir_10m"
    CANONICAL_UNITS['wind_from_direction_10m'] = "degree" # from North=0 / East=90"










