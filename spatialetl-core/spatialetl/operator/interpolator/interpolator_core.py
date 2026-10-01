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

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, UTC

import numpy as np
from numpy import int8, int16, int32, int64
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, CloughTocher2DInterpolator
from scipy.interpolate import interp1d
from scipy.spatial import qhull

from spatialetl.utils.logger import logging
from spatialetl.utils.timing import timing


def interp_2d_weights(gridX, gridY, map, interpolator="scipy"):
    """
    Interpolate the 2d weights of the source grid

    :param gridX:
    :param gridY:
    :param current_thread (int): Current thread
    :return current_thread, Delaunay triangulation
    """

    logging.debug(f"[InterpolatorCore][horizontal_interpolation()] computing weights with '{interpolator}'")

    gridX = np.ma.filled(gridX, fill_value=-9999.)
    gridY = np.ma.filled(gridY, fill_value=-9999.)

    if interpolator == "cgal":
        try:
            from spatialetl.operator.interpolator.interpcgal.nninterpol import triangulate

            threads_count = np.shape(map)[0]
            tri = np.empty([threads_count], dtype=object)

            with ThreadPoolExecutor(max_workers=threads_count) as executor:
                futures = []
                for current_thread in range(0, executor._max_workers):
                    if gridX.ndim == 1 and gridY.ndim == 1:
                        local_grid_x, local_grid_y = np.meshgrid(gridX[map[current_thread]["src_global_x_overlap"]],
                                                                 gridY[map[current_thread]["src_global_y_overlap"]])
                    else:
                        local_grid_x = gridX[map[current_thread]["src_global_y_overlap"],
                        map[current_thread]["src_global_x_overlap"]]
                        local_grid_y = gridY[map[current_thread]["src_global_y_overlap"],
                        map[current_thread]["src_global_x_overlap"]]

                    points = np.array([local_grid_x.flatten(), local_grid_y.flatten()]).T

                    futures.append(
                        executor.submit(triangulate, points))

                for current_thread, f in enumerate(futures):
                    tri[current_thread] = f.result()

            return tri
        except ModuleNotFoundError:
            logging.warning(f"Unable to find the cgal interpolator, we switch to scipy")
            interpolator = "scipy"

    if interpolator == "scipy":
        if gridX.ndim == 1 and gridY.ndim == 1:
            local_grid_x, local_grid_y = np.meshgrid(gridX, gridY)
        else:
            local_grid_x = gridX
            local_grid_y = gridY

        points = np.array([local_grid_x.flatten(), local_grid_y.flatten()]).T
        return qhull.Delaunay(points)
    else:
        raise ValueError(f"Unable to find interpolator {interpolator}")


def resample_2d_to_grid(tri, newX, newY, data, method, map, interpolator="scipy",is_binary_data=False):
    """
    2D resampling function
    """
    logging.debug(f"[InterpolatorCore][horizontal_interpolation()] Starting interpolation with '{interpolator}'")

    if data.dtype == int8 or data.dtype == int16 or data.dtype == int32 or data.dtype == int64:
        fill_value = -9999
    else:
        fill_value = 9.96921e+36

    if interpolator == "cgal":
        try:
            from spatialetl.operator.interpolator.interpcgal.nninterpol import nninterpol
            threads_count = np.shape(map)[0]
            local_data = np.zeros([np.shape(newY)[0], np.shape(newX)[0]])
            local_data[:] = fill_value

            if method == "linear":
                pass
            elif method == "nearest" or method == "cubic" :
                raise ValueError(
                    f"Interpolation method '{method}' is not implemented in cgal interpolator'. Use 'linear' method")
            else:
                raise ValueError(
                    f"Unable to decode interpolation method '{method}'. Try 'linear' or 'nearest' or 'cubic'")

            with ThreadPoolExecutor(max_workers=threads_count) as executor:
                futures = []
                for current_thread in range(0, executor._max_workers):
                    dst_local_grid_x, dst_local_grid_y = np.meshgrid(newX[map[current_thread]["dst_parent_x_overlap"]],
                                                                     newY[map[current_thread]["dst_parent_y_overlap"]])

                    values = data[map[current_thread]["src_global_y_overlap"],
                    map[current_thread]["src_global_x_overlap"]].flatten().T

                    if is_binary_data:
                        # To avoid doubtful linear interpolation with binary data,
                        # we change enlarge the boundaries to -9999 / 9999
                        values[values == 0] = -9999
                        values[values == 1] = 9999

                    futures.append(
                        executor.submit(nninterpol, tri[current_thread], values, dst_local_grid_x, dst_local_grid_y,
                                        fill_value))

                for current_thread, f in enumerate(futures):
                    data = f.result()
                    slice_data = np.array(data)
                    local_data[map[current_thread]["dst_parent_y"], map[current_thread]["dst_parent_x"]] = slice_data[
                        map[current_thread]["dst_local_y"],
                        map[current_thread]["dst_local_x"]]

            # We force the result's type to be the same as the original one
            local_data = local_data.astype(values.dtype)

            if is_binary_data:
                # We avoid doubtful linear interpolation with binary data,
                # we compute binary data from 0
                local_data[local_data <= 0] = 0
                local_data[local_data > 0] = 1

            return local_data
        except ModuleNotFoundError:
            logging.warning(f"Unable to find the cgal interpolator, we switch to scipy")
            interpolator = "scipy"

    if interpolator == "scipy":
        # Target grid
        dst_local_grid_x, dst_local_grid_y = np.meshgrid(newX, newY)

        # Values
        values = data.flatten()

        if method == "linear":

            if is_binary_data:
                # To avoid doubtful linear interpolation with binary data,
                # we change enlarge the boundaries to -9999 / 9999
                values[values == 0] = -9999
                values[values == 1] = 9999

            ip = LinearNDInterpolator(tri, values, fill_value=fill_value,
                                      rescale=False)
        elif method == "nearest":
            ip = NearestNDInterpolator(tri, values, rescale=False)
        elif method == "cubic":

            if is_binary_data:
                # To avoid doubtful linear interpolation with binary data,
                # we change enlarge the boundaries to -9999 / 9999
                values[values == 0] = -9999
                values[values == 1] = 9999

            ip = CloughTocher2DInterpolator(tri, values, rescale=False)
        else:
            raise ValueError(f"Unable to decode interpolation method '{method}'. Try 'linear' or 'nearest' or 'cubic'")

        data = ip((dst_local_grid_x, dst_local_grid_y))

        # We force the result's type to be the same as the original one
        data = data.astype(values.dtype)

        if is_binary_data and (method == "linear" or method == "cubic"):
            # We avoid doubtful linear interpolation with binary data,
            # we compute binary data from 0
            data[data <= 0] = 0
            data[data > 0] = 1

        return data
    else:
        raise ValueError(f"Unable to find interpolator {interpolator}")


def vertical_interpolation(sourceAxis, targetAxis, data, method, extrapolate=False):
    logging.debug("[InterpolatorCore][vertical_interpolation()] Looking for water depth : " + str(
       targetAxis[0]) + " m with method '" + str(method) + "'.")
    logging.debug("[InterpolatorCore][vertical_interpolation()] Source Axis contains: " + str(sourceAxis))
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Candidates values are: " + str(data))
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Target Axis contains: " + str(targetAxis))
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Method: " + str(method))
    #logging.debug("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")

    if method == "mean":
        return np.mean(data)

    elif method == "nearest":
        array = np.asarray(sourceAxis)
        nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
        return data[nearest_index_t]
    elif method == "linear" or method == "cubic":
        try:
            f = interp1d(sourceAxis, data, kind=method, bounds_error=True)
            return f(targetAxis)
        except ValueError as ex:
            logging.warning("[InterpolatorCore][vertical_interpolation()] Error: " + str(ex))
            logging.warning(
                "[InterpolatorCore][vertical_interpolation()] This error may occur when you ask for a water depth out of range: " + str(
                    targetAxis[0]) + " m")
            logging.warning(
                "[InterpolatorCore][vertical_interpolation()] We found these water depth candidates: " + str(
                    sourceAxis))
            logging.warning(
                "[InterpolatorCore][vertical_interpolation()] To avoid this error, you can change your zbox range or use another vertical interpolation method.")
            if extrapolate:
                logging.warning(
                    "[InterpolatorCore][vertical_interpolation()] We continue by using an extrapolation method.")
                logging.warning("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")
                f = interp1d(sourceAxis, data, kind=method, fill_value="extrapolate")
                return f(targetAxis)
            else:
                logging.warning("[InterpolatorCore][vertical_interpolation()] We continue by using the nearest method.")
                logging.warning("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")
                array = np.asarray(sourceAxis)
                nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
                return data[nearest_index_t]
    else:
        raise ValueError("Unable to decode vertical interpolation method : " + str(method))


def temporal_1d_interpolation(sourceAxis, targetAxis, data, method, extrapolate=False):
    logging.debug("[InterpolatorCore][temporal_interpolation()] Looking for time : " + str(
        datetime.utcfromtimestamp(targetAxis[0])) + " with method '" + str(method) + "'.")
    for time in sourceAxis:
        logging.debug(
            "[InterpolatorCore][temporal_interpolation()] Source Axis contains: " + str(
                datetime.utcfromtimestamp(time)))

    if method is None:
        return np.nan

    elif method == "mean":
        return np.mean(data)

    elif method == "nearest":

        array = np.asarray(sourceAxis)
        nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
        return data[nearest_index_t]

    elif method == "linear" or method == "cubic":

        try:
            f = interp1d(sourceAxis, data, kind=method, bounds_error=True)
            return f(targetAxis)
        except ValueError as ex:
            logging.warning("[InterpolatorCore][temporal_interpolation()] Error: " + str(ex))
            logging.warning(
                "[InterpolatorCore][temporal_interpolation()] This error may occur when you ask for a datetime  out of range : " + str(
                    datetime.utcfromtimestamp(targetAxis[0])))
            logging.warning(
                "[InterpolatorCore][temporal_interpolation()] To avoid this error, you can change your time range or use another time interpolation method.")
            if extrapolate:
                logging.warning(
                    "[InterpolatorCore][temporal_interpolation()] We continue by using an extrapolation method.")
                logging.warning("[InterpolatorCore][temporal_interpolation()] ----------------------------------------")
                f = interp1d(sourceAxis, data, kind=method, fill_value="extrapolate")
                return f(targetAxis)
            else:
                logging.warning("[InterpolatorCore][temporal_interpolation()] We continue by using the nearest method.")
                logging.warning("[InterpolatorCore][temporal_interpolation()] ----------------------------------------")
                array = np.asarray(sourceAxis)
                nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
                return data[nearest_index_t]
    else:
        raise ValueError("Unable to decode temporal interpolation method : " + str(method))


def temporal_2d_interpolation(source_axis, target_time, data,method,time_delta):
    logging.debug(f"[InterpolatorCore][temporal_interpolation()] Starting interpolation with method '{method}'.")

    if method == "nearest":
        nearest_t_index = (np.abs(source_axis - target_time)).argmin()
        if target_time - source_axis[nearest_t_index] == 0.0 or abs(
                target_time - source_axis[nearest_t_index]) < time_delta:

            logging.debug("[TimeCoverage][find_time_index()] Nearest datetime found : " + str(
                source_axis[nearest_t_index]))

            return data[nearest_t_index,:]

        else:
            raise ValueError("'" + str(
                                          target_time) + "' not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                                          time_delta) + ") is too small or the date is out the range.")
    elif method == "mean":
        return np.mean(data, axis=0)
    else:
        raise ValueError("Unable to decode temporal interpolation method : " + str(method))
