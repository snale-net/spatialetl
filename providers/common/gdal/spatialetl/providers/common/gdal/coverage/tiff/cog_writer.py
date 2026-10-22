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
from osgeo import gdal
from osgeo import osr

from spatialetl.coverage import LevelCoverage
from spatialetl.coverage.time_coverage import TimeCoverage
from spatialetl.coverage.time_level_coverage import TimeLevelCoverage
from spatialetl.exception.coverage_error import CoverageError
from spatialetl.providers.common.gdal.coverage.tiff.default_writer import DefaultWriter as TiffWriter
from spatialetl.utils.logger import logging


class COGWriter(TiffWriter):
    """
    Writer to write Cloud Optimized GeoTiff (COG)
    """

    def __init__(self, cov, myFile):
        TiffWriter.__init__(self, cov, myFile);
        self.driver = gdal.GetDriverByName("MEM")

    def _create_file(self,
                     variable_name,
                     time=None,
                     depth=None,
                     variable_type=gdal.GDT_Float32,
                     variable_fill_value=9.96921e+36):
        if self.coverage.rank == 0:
            if time and depth:
                file = self.driver.Create("",
                                          int(self.rows), int(self.cols), 1, variable_type)
            elif depth:
                file = self.driver.Create("",
                                          int(self.rows), int(self.cols), 1, variable_type)
            elif time:
                file = self.driver.Create("",
                                          int(self.rows), int(self.cols), 1, variable_type)
            else:
                file = self.driver.Create("",
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
            logging.info(f"[COGWriter] Writing variable '{variable_long_name}'")
        fn = getattr(self.coverage, function_name)

        self._gathering_var(fn(), var)

        if self.coverage.rank == 0:
            file.GetRasterBand(1).WriteArray(var)

            self._write_cog(
                file,
                variable_name=variable_name,
            )

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
                        f"[COGWriter] Writing variable '{variable_long_name}' at time '{time}'")

                file = self._create_file(
                    variable_name,
                    time=time,
                    variable_type=variable_type,
                    variable_fill_value=variable_fill_value)

                self._gathering_var(fn(time), var, time_index=time_index)

                if self.coverage.rank == 0:
                    file.GetRasterBand(1).WriteArray(var[time_index])

                    self._write_cog(
                        file,
                        variable_name=variable_name,
                        time=time)

                time_index += 1
        else:
            raise CoverageError("COGWriter",
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
                        f"[COGWriter] Writing variable '{variable_name}' at time '{time}'")

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

                if self.coverage.rank == 0:
                    ufile.GetRasterBand(1).WriteArray(ucomp[time_index])
                    vfile.GetRasterBand(1).WriteArray(vcomp[time_index])

                    self._write_cog(
                        ufile,
                        variable_name=eastward_variable_name,
                        time=time)

                    self._write_cog(
                        vfile,
                        variable_name= northward_variable_name,
                        time=time)

                time_index += 1
        else:
            raise CoverageError("COGWriter",
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
                    f"[COGWriter] Writing variable '{variable_long_name}'")

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

                if self.coverage.rank == 0:
                    file.GetRasterBand(1).WriteArray(var[level_index])

                    self._write_cog(
                        file,
                        variable_name,
                        depth=level)

                level_index += 1
        else:
            raise CoverageError("COGWriter",
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
                        f"[COGWriter] Writing variable '{variable_long_name}' at time '{time}'")
                level_index = 0
                for level in self.coverage.read_axis_z():
                    file = self._create_file(
                        variable_name,
                        time=time,
                        depth=level,
                        variable_type=variable_type,
                        variable_fill_value=variable_fill_value)

                    self._gathering_var(fn(time, level), var, time_index=time_index, level_index=level_index)

                    if self.coverage.rank == 0:
                        file.GetRasterBand(1).WriteArray(var[time_index, level_index])

                        self._write_cog(
                            file,
                            variable_name,
                            time=time,
                            depth=level)

                    level_index += 1

                time_index += 1
        else:
            raise CoverageError("COGWriter", "The given coverage is not an instance of 'TimeLevelCoverage'")

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
                        f"[COGWriter] Writing variable '{variable_name}' at time '{time}'")

                level_index = 0
                for level in self.coverage.read_axis_z():
                    local_data_u, local_data_v = fn(time, level)

                    ufile = self._create_file(
                        eastward_variable_name,
                        time=time,
                        depth=level,
                        variable_type=variable_type,
                        variable_fill_value=variable_fill_value)

                    vfile = self._create_file(
                        northward_variable_name,
                        time=time,
                        depth=level,
                        variable_type=variable_type,
                        variable_fill_value=variable_fill_value)

                    self._gathering_var(local_data_u, ucomp, time_index=time_index, level_index=level_index)
                    self._gathering_var(local_data_v, vcomp, time_index=time_index, level_index=level_index)

                    if self.coverage.rank == 0:
                        ufile.GetRasterBand(1).WriteArray(ucomp[time_index, level_index])
                        vfile.GetRasterBand(1).WriteArray(vcomp[time_index, level_index])

                        self._write_cog(
                            ufile,
                            variable_name=eastward_variable_name,
                            time=time,
                            depth=level)

                        self._write_cog(
                            vfile,
                            variable_name=northward_variable_name,
                            time=time,
                            depth=level)

                    level_index += 1
                time_index += 1
        else:
            raise CoverageError("COGWriter",
                                "The given coverage is not an instance of 'TimeLevelCoverage'")

    def _write_cog(self,
                   ds,
                   variable_name,
                   time=None,
                   depth=None):
        if self.coverage.rank == 0:
            if time and depth:
                file = os.path.join(
                    self.filename,
                    f"{time.strftime('%Y%m%d_%H%M%S')}_{depth}m_{variable_name}.tiff"
                )
            elif depth:
                file = os.path.join(
                    self.filename,
                    f"{depth}m_{variable_name}.tiff"
                )
            elif time:
                file = os.path.join(
                    self.filename,
                    f"{time.strftime('%Y%m%d_%H%M%S')}_{variable_name}.tiff"
                )
            else:
                file = os.path.join(
                    self.filename,
                    f"{variable_name}.tiff"
                )

            # Write COG file
            gdal.Translate(
                file,
                ds,
                format="COG",
                creationOptions=[
                    "COMPRESS=DEFLATE",
                    "PREDICTOR=FLOATING_POINT",
                    "BLOCKSIZE=512",
                    "OVERVIEWS=AUTO",
                    "BIGTIFF=IF_SAFER",
                ],
            )

        ds = None
