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

from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32, float64, int32

from spatialetl.point.io.multi_point_writer import MultiPointWriter
from spatialetl.point.level_multi_point import LevelMultiPoint
from spatialetl.point.time_level_multi_point import TimeLevelMultiPoint
from spatialetl.point.time_multi_point import TimeMultiPoint
from spatialetl.utils.logger import logging
from spatialetl.utils.variable_definition import VariableDefinition


class DefaultWriter(MultiPointWriter):

    def __init__(self, p, myFile):
        MultiPointWriter.__init__(self, p, myFile)
        format = 'NETCDF4_CLASSIC'

        self.ncfile = Dataset(self.filename, 'w', format=format)
        self.ncfile.description = 'Generated with spatialETL'

        self.ncfile.data_source = str(self.points.data_source)
        self.ncfile.meta_data = str(self.points.meta_data)

        # Geo-points dimension
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['point'], self.points.get_nb_points())
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['point'], int32,
                                         (VariableDefinition.VARIABLE_NAME['point'],))
        var.long_name = VariableDefinition.LONG_NAME['point']
        var.standard_name = VariableDefinition.STANDARD_NAME['point']
        var.axis = "X";
        var.units = VariableDefinition.CANONICAL_UNITS['point'];
        var[:] = range(0, self.points.get_nb_points());

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['latitude'], float32,
                                         (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['latitude']
        var.standard_name = VariableDefinition.STANDARD_NAME['latitude']
        var.valid_min = -90.0;
        var.valid_max = 90.0;
        var.units = VariableDefinition.CANONICAL_UNITS['latitude']
        var[:] = self.points.read_axis_y()

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['longitude'], float32,
                                         (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['longitude']
        var.standard_name = VariableDefinition.STANDARD_NAME['longitude']
        var.valid_min = -180.0;
        var.valid_max = 180.0;
        var.units = VariableDefinition.CANONICAL_UNITS['longitude']
        var[:] = self.points.read_axis_x()

        if (isinstance(self.points, TimeMultiPoint) or isinstance(self.points, TimeLevelMultiPoint)):
            # Time dimension
            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['time'],
                                        self.points.get_t_size(type="target_global"))
            times = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['time'], float64,
                                               (VariableDefinition.VARIABLE_NAME['time'],))
            times.units = 'seconds since 1970-01-01 00:00:00'
            times.calendar = 'gregorian'
            times.standard_name = 'time'
            times.axis = 'T'
            times.conventions = "UTC time"

            times[:] = date2num(self.points.read_axis_t(),
                                units=times.units,
                                calendar=times.calendar)

        if (isinstance(self.points, LevelMultiPoint) or isinstance(self.points, TimeLevelMultiPoint)):

            if self.points.is_sigma_coordinate():
                raise ValueError("This writer supports only regular vertical axis.")

            # Depth dimension
            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['depth'], self.points.get_z_size())
            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['depth'], float32,
                                             (VariableDefinition.VARIABLE_NAME['depth'],))
            var.standard_name = VariableDefinition.STANDARD_NAME['depth']
            var.long_name = VariableDefinition.LONG_NAME['depth']
            var.positive = "down";
            var.axis = 'Z'
            var.units = VariableDefinition.CANONICAL_UNITS['depth'];
            var[:] = self.points.read_axis_z()

    def close(self):
        self.ncfile.close()

    def write_variable_time(self):
        times = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['time'], float64,
                                           (VariableDefinition.VARIABLE_NAME['point'],))
        times.units = 'seconds since 1970-01-01 00:00:00'
        times.calendar = 'gregorian'
        times.standard_name = 'time'
        times.conventions = "UTC time"

        times[:] = date2num(self.points.read_variable_time(), units=times.units, calendar=times.calendar)

    #################
    # HYDRO
    # 2D
    #################

    def write_variable_bathymetry(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['bathymetry'], float32,
                                         (VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['bathymetry']
        var.standard_name = VariableDefinition.STANDARD_NAME['bathymetry']
        var.units = VariableDefinition.CANONICAL_UNITS['bathymetry']
        var.positive = "down";

        logging.info(
            '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['bathymetry']) + '\'')
        var[:] = self.points.read_variable_bathymetry()

    def write_variable_sea_surface_height_above_mean_sea_level(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_height_above_mean_sea_level'],
                                         float32, (
                                             VariableDefinition.VARIABLE_NAME['time'],
                                             VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_height_above_mean_sea_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_mean_sea_level']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']) + '\' at time \'' + str(
                    time) + '\'')

            var[
                time_index:time_index + 1, :] = self.points.read_variable_sea_surface_height_above_mean_sea_level_at_time(
                time)
            time_index += 1

    def write_variable_water_volume_transport_into_sea_water_from_rivers(self):
        var = self.ncfile.createVariable(
            VariableDefinition.VARIABLE_NAME['water_volume_transport_into_sea_water_from_rivers'], float32, (
                VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']),
            fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['water_volume_transport_into_sea_water_from_rivers']
        var.standard_name = VariableDefinition.STANDARD_NAME['water_volume_transport_into_sea_water_from_rivers']
        var.units = VariableDefinition.CANONICAL_UNITS['water_volume_transport_into_sea_water_from_rivers']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['water_volume_transport_into_sea_water_from_rivers']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME[
                        'water_volume_transport_into_sea_water_from_rivers']) + '\' at time \'' + str(
                    time) + '\'')

            var[
                time_index:time_index + 1, :] = self.points.read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(
                time)
            time_index += 1

    def write_variable_barotropic_sea_water_speed(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_sea_water_speed'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['barotropic_sea_water_speed']
        var.standard_name = VariableDefinition.STANDARD_NAME['barotropic_sea_water_speed']
        var.units = VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_speed']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['barotropic_sea_water_speed']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['barotropic_sea_water_speed']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_speed_at_time(time)
            time_index += 1

    def write_variable_barotropic_sea_water_from_direction(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_sea_water_from_direction'],
                                         float32, (VariableDefinition.VARIABLE_NAME['time'],
                                                   VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction']
        var.standard_name = VariableDefinition.STANDARD_NAME['barotropic_sea_water_from_direction']
        var.units = VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_from_direction']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_from_direction_at_time(
                time)
            time_index += 1

    def write_variable_barotropic_sea_water_to_direction(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_sea_water_to_direction'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'],
                                          VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction']
        var.standard_name = VariableDefinition.STANDARD_NAME['barotropic_sea_water_to_direction']
        var.units = VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_to_direction']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_to_direction_at_time(
                time)
            time_index += 1

    def write_variable_ocean_tracer_residence_time(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['ocean_tracer_residence_time'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'],
                                          VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['ocean_tracer_residence_time']
        var.standard_name = VariableDefinition.STANDARD_NAME['ocean_tracer_residence_time']
        var.units = VariableDefinition.CANONICAL_UNITS['ocean_tracer_residence_time']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['ocean_tracer_residence_time']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['ocean_tracer_residence_time']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_ocean_tracer_residence_time_at_time(
                time)
            time_index += 1

    #################
    # HYDRO
    # Sea Surface
    #################

    def write_variable_sea_surface_temperature(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_temperature'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_temperature'];
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_temperature']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_temperature']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_surface_temperature_at_time(time)
            time_index += 1

    def write_variable_sea_surface_salinity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_salinity'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_salinity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_salinity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_salinity']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_surface_salinity_at_time(time)
            time_index += 1

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_pressure_at_sea_water_surface'],
                                         float32, (
                                             VariableDefinition.VARIABLE_NAME['time'],
                                             VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_pressure_at_sea_water_surface']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_pressure_at_sea_water_surface']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1,
            :] = self.points.read_variable_sea_water_pressure_at_sea_water_surface_at_time(time)
            time_index += 1

    def write_variable_sea_water_speed_at_sea_water_surface(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_speed_at_sea_water_surface'],
                                         float32, (
                                             VariableDefinition.VARIABLE_NAME['time'],
                                             VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_speed_at_sea_water_surface']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_speed_at_sea_water_surface']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_speed_at_sea_water_surface']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_speed_at_sea_water_surface']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_speed_at_sea_water_surface']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_water_speed_at_sea_water_surface_at_time(
                time)
            time_index += 1

    def write_variable_sea_water_from_direction_at_sea_water_surface(self):
        var = self.ncfile.createVariable(
            VariableDefinition.VARIABLE_NAME['sea_water_from_direction_at_sea_water_surface'],
            float32, (VariableDefinition.VARIABLE_NAME['time'],
                      VariableDefinition.VARIABLE_NAME['point']),
            fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_from_direction_at_sea_water_surface']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_from_direction_at_sea_water_surface']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_from_direction_at_sea_water_surface']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_from_direction_at_sea_water_surface']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME[
                        'sea_water_from_direction_at_sea_water_surface']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1,
            :] = self.points.read_variable_sea_water_from_direction_at_sea_water_surface_at_time(
                time)
            time_index += 1

    def write_variable_sea_water_to_direction_at_sea_water_surface(self):
        var = self.ncfile.createVariable(
            VariableDefinition.VARIABLE_NAME['sea_water_to_direction_at_sea_water_surface'],
            float32,
            (VariableDefinition.VARIABLE_NAME['time'],
             VariableDefinition.VARIABLE_NAME['point']),
            fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_to_direction_at_sea_water_surface']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_to_direction_at_sea_water_surface']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_to_direction_at_sea_water_surface']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_to_direction_at_sea_water_surface']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME[
                        'sea_water_to_direction_at_sea_water_surface']) + '\' at time \'' + str(
                    time) + '\'')

            var[
                time_index:time_index + 1, :] = self.points.read_variable_sea_water_to_direction_at_sea_water_surface_at_time(
                time)
            time_index += 1

    #################
    # HYDRO
    # Ground level
    #################

    def write_variable_sea_water_temperature_at_ground_level(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_temperature_at_ground_level'],
                                         float32, (
                                             VariableDefinition.VARIABLE_NAME['time'],
                                             VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level'];
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_temperature_at_ground_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_temperature_at_ground_level']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_water_temperature_at_ground_level_at_time(
                time)
            time_index += 1

    def write_variable_sea_water_salinity_at_ground_level(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_salinity_at_ground_level'],
                                         float32, (
                                             VariableDefinition.VARIABLE_NAME['time'],
                                             VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_salinity_at_ground_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_salinity_at_ground_level']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_water_salinity_at_ground_level_at_time(
                time)
            time_index += 1

    def write_variable_sea_water_speed_at_ground_level(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_speed_at_ground_level'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_speed_at_ground_level']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_speed_at_ground_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_speed_at_ground_level']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_speed_at_ground_level']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_speed_at_ground_level']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_water_speed_at_ground_level_at_time(time)
            time_index += 1

    def write_variable_sea_water_from_direction_at_ground_level(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_from_direction_at_ground_level'],
                                         float32, (VariableDefinition.VARIABLE_NAME['time'],
                                                   VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_from_direction_at_ground_level']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_from_direction_at_ground_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_from_direction_at_ground_level']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_from_direction_at_ground_level']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_from_direction_at_ground_level']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1,
            :] = self.points.read_variable_sea_water_from_direction_at_ground_level_at_time(
                time)
            time_index += 1

    def write_variable_sea_water_to_direction_at_ground_level(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_to_direction_at_ground_level'],
                                         float32,
                                         (VariableDefinition.VARIABLE_NAME['time'],
                                          VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_to_direction_at_ground_level']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_to_direction_at_ground_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_to_direction_at_ground_level']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_to_direction_at_ground_level']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_to_direction_at_ground_level']) + '\' at time \'' + str(
                    time) + '\'')

            var[
                time_index:time_index + 1, :] = self.points.read_variable_sea_water_to_direction_at_ground_level_at_time(
                time)
            time_index += 1

    #################
    # HYDRO
    # 3D
    #################

    def write_variable_sea_water_temperature(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_temperature'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'],
            VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_temperature']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_temperature']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_temperature']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_temperature']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info('[DefaultWriter] Writing variable \'' + str(
                VariableDefinition.LONG_NAME['sea_water_temperature']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_temperature_at_time_and_depth(time, depth)
                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_salinity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_salinity'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'],
            VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_salinity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_salinity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_salinity']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_salinity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_salinity']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_salinity_at_time_and_depth(time, depth)
                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_density(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_density'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'],
            VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_density']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_density']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_density']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_density']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_density']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_density_at_time_and_depth(time, depth)
                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_electrical_conductivity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_electrical_conductivity'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'],
                                          VariableDefinition.VARIABLE_NAME['depth'],
                                          VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_electrical_conductivity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_electrical_conductivity']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']) + '\' at time \'' + str(
                    time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_conductivity_at_time_and_depth(time, depth)
                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_turbidity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_turbidity'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'],
            VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_turbidity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_turbidity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_turbidity']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_turbidity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_turbidity']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_turbidity_at_time_and_depth(time, depth)
                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_baroclinic_sea_water_velocity(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'],
                                          float32, (VariableDefinition.VARIABLE_NAME['time'],
                                                    VariableDefinition.VARIABLE_NAME['depth'],
                                                    VariableDefinition.VARIABLE_NAME['point'],),
                                          fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['baroclinic_eastward_sea_water_velocity']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_eastward_sea_water_velocity']
        ucur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity']
        ucur.comment = "cur=sqrt(U**2+V**2)";

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_northward_sea_water_velocity'],
                                          float32, (VariableDefinition.VARIABLE_NAME['time'],
                                                    VariableDefinition.VARIABLE_NAME['depth'],
                                                    VariableDefinition.VARIABLE_NAME['point']),
                                          fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['baroclinic_northward_sea_water_velocity']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_northward_sea_water_velocity']
        vcur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_northward_sea_water_velocity']
        vcur.comment = "cur=sqrt(U**2+V**2)";

        logging.info('[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\'\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\' at time \'' + str(
                    time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                cur = self.points.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(time, depth)
                ucur[time_index:time_index + 1, z_index:z_index + 1] = cur[0]
                vcur[time_index:time_index + 1, z_index:z_index + 1] = cur[1]
                z_index = z_index + 1

            time_index += 1

    #################
    # WAVES
    # Sea Surface
    #################

    #################
    # WAVES
    # Ground level
    #################

    #################
    # WAVES
    # Momentum flux
    #################

    #################
    # METEO
    # 2D
    #################

    #################
    # METEO
    # Surface air
    #################

    def write_variable_rainfall_amount(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['rainfall_amount'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['rainfall_amount']
        var.standard_name = VariableDefinition.STANDARD_NAME['rainfall_amount']
        var.units = VariableDefinition.CANONICAL_UNITS['rainfall_amount']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['rainfall_amount']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['rainfall_amount']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_rainfall_amount_at_time(time)
            time_index += 1

    def write_variable_surface_air_pressure(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['surface_air_pressure'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['surface_air_pressure']
        var.standard_name = VariableDefinition.STANDARD_NAME['surface_air_pressure']
        var.units = VariableDefinition.CANONICAL_UNITS['surface_air_pressure']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['surface_air_pressure']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_air_pressure']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_surface_air_pressure_at_time(time)
            time_index += 1

    #################
    # METEO
    # At 10 m
    #################

    def write_variable_wind_speed_10m(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_speed_10m'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['wind_speed_10m']
        var.standard_name = VariableDefinition.STANDARD_NAME['wind_speed_10m']
        var.units = VariableDefinition.CANONICAL_UNITS['wind_speed_10m']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['wind_speed_10m']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_speed_10m']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_speed_10m_at_time(time)
            time_index += 1

    def write_variable_wind_from_direction_10m(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_from_direction_10m'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['wind_from_direction_10m']
        var.standard_name = VariableDefinition.STANDARD_NAME['wind_from_direction_10m']
        var.units = VariableDefinition.CANONICAL_UNITS['wind_from_direction_10m']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_from_direction_10m_at_time(time)
            time_index += 1

    def write_variable_wind_to_direction_10m(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_to_direction_10m'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['wind_to_direction_10m']
        var.standard_name = VariableDefinition.STANDARD_NAME['wind_to_direction_10m']
        var.units = VariableDefinition.CANONICAL_UNITS['wind_to_direction_10m']

        logging.info('[DefaultWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['wind_to_direction_10m']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.info(
                '[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_to_direction_10m']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_to_direction_10m_at_time(time)
            time_index += 1
