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
import inspect
import os
from importlib.util import find_spec

import numpy as np
from array_split import shape_split

from spatialetl.exception.not_found_in_rank_error import NotFoundInRankError
from spatialetl.operator.interpolator.interpolator_core import resample_2d_to_grid, interp_2d_weights
from spatialetl.utils.distance import distance_on_unit_sphere
from spatialetl.utils.logger import logging

mpi_lib = find_spec("mpi4py")
MPI_FOUND = mpi_lib is not None
if MPI_FOUND:
    from mpi4py import MPI


class Coverage(object):
    """
    The Coverage class represents a spatial coverage on the horizontal plane.

    The points representing this coverage can be aligned on a regular grid (x,y) or on an irregular grid ((x1,y1),(x2,y2)).
    Depending on the grid type, the axis reading functions will return one or two-dimensional arrays.
    To avoid loading the entire file into memory, the coverage contains a pointer to a reader.
    Layers are therefore read on demand from the file.

    Note: The axes are always reversed in the arrays because of NetCDF.
    So the y-axis first, then the x-axis. Example: [y,x]

    Examples
    --------
    >>> coverage = Coverage(myReader, bbox=[-5, 5, -5, 5], resolution_x=0.1, resolution_y=0.1)
    """

    HORIZONTAL_INTERPOLATOR = "scipy"
    HORIZONTAL_INTERPOLATION_METHOD = "linear"
    HORIZONTAL_OVERLAPING_SIZE = 2

    def __init__(self, reader, bbox=None, resolution_x=None, resolution_y=None, parallel=True,
                 nb_thread: int = os.cpu_count() - 1):
        """
        Initialize the Coverage object with a file reader and optional bounding box and resolution.
        
        Parameters
        ----------
        
        reader : object
            File reader instance.
        bbox : list, optional
            Bounding box coordinates [xmin, xmax, ymin, ymax].
        resolution_x : float, optional
            Resolution along the x-axis.
        resolution_y : float, optional
            Resolution along the y-axis.

        Examples
        --------
        >>> coverage = Coverage(reader, bbox=[0, 10, 0, 10], resolution_x=1.0, resolution_y=1.0)
        """

        self.reader = reader;

        if parallel:
            if MPI_FOUND:
                # Parallel MPI
                self.comm = MPI.COMM_WORLD
                self.size = self.comm.Get_size()
                self.rank = self.comm.Get_rank()
            else:
                # No parallel MPI
                self.comm = None
                self.parallel_map = None
                self.size = 1
                self.rank = 0

            # Parallel multithreading
            self.threads_number = nb_thread
        else:
            # No parallel MPI
            self.comm = None
            self.parallel_map = None
            self.size = 1
            self.rank = 0

            # No parallel multithreading
            self.threads_number = 1

        self.parallel_map = np.empty(self.size, dtype=object)

        self.source_regular_grid = self.reader.is_regular_grid()
        self.target_regular_grid = self.source_regular_grid
        self.horizontal_resampling = False

        self.source_global_x_size = self.reader.get_x_size()
        self.source_global_y_size = self.reader.get_y_size()
        self.source_global_axis_x = self.reader.read_axis_x(0, self.source_global_x_size, 0,
                                                            self.source_global_y_size)
        self.source_global_axis_y = self.reader.read_axis_y(0, self.source_global_x_size, 0,
                                                            self.source_global_y_size)
        self.source_global_tri = None

        # We adjust the computation based on the bbox
        if bbox is None:
            # we compute the destination grid
            Ymin = np.min(self.source_global_axis_y)
            Ymax = np.max(self.source_global_axis_y)
            Xmin = np.min(self.source_global_axis_x)
            Xmax = np.max(self.source_global_axis_x)
        else:
            if self.check_bbox_validity(bbox) is False:
                raise ValueError("Your Bbox is not valid (e.g [Xmin,Xmax,Ymin,Ymax] or is outside the coverage")
            Ymin = bbox[2]
            Ymax = bbox[3]
            Xmin = bbox[0]
            Xmax = bbox[1]

        if self.is_regular_grid(type="source"):

            idx = np.where((self.source_global_axis_x >= Xmin) &
                           (self.source_global_axis_x <= Xmax))

            xmin = np.min(idx[0])
            xmax = np.max(idx[0]) + 1

            idx = np.where((self.source_global_axis_y >= Ymin) &
                           (self.source_global_axis_y <= Ymax))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1

            self.target_global_axis_x = self.source_global_axis_x[xmin:xmax]
            self.target_global_x_size = xmax - xmin
            self.target_global_axis_y = self.source_global_axis_y[ymin:ymax]
            self.target_global_y_size = ymax - ymin

        else:

            idx = np.where((self.source_global_axis_x >= Xmin) &
                           (self.source_global_axis_x <= Xmax) &
                           (self.source_global_axis_y >= Ymin) &
                           (self.source_global_axis_y <= Ymax))

            if len(idx[0]) == 0:
                raise ValueError("No values found")

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1
            xmin = np.min(idx[1])
            xmax = np.max(idx[1]) + 1

            self.target_global_axis_x = self.source_global_axis_x[ymin:ymax, xmin:xmax]
            self.target_global_x_size = xmax - xmin
            self.target_global_axis_y = self.source_global_axis_y[ymin:ymax, xmin:xmax]
            self.target_global_y_size = ymax - ymin

        # Compute the destination grid
        self.target_global_res_x = None
        self.target_global_res_y = None

        if resolution_x is not None and resolution_y is not None:

            self.horizontal_resampling = True
            self.target_regular_grid = True

            res = np.mean([resolution_x, resolution_y])
            self.target_global_res_x = res
            self.target_global_res_y = res

            self.target_global_axis_x = np.arange(Xmin, Xmax, res)
            self.target_global_axis_y = np.arange(Ymin, Ymax, res)

            self.target_global_x_size = len(self.target_global_axis_x)
            self.target_global_y_size = len(self.target_global_axis_y)

            if type(self) == Coverage and self.rank == 0:
                logging.info(
                    '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                        self.source_global_y_size) + ")")
                logging.info(
                    '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                        self.target_global_y_size) + ")")

        if type(self) == Coverage:
            self.__init_parallel_map()

            if self.rank == 0 and MPI_FOUND:
                logging.debug("MPI map:")

            if MPI_FOUND:
                logging.debug(f"{'-' * 10} MPI rank n° {self.rank} {'-' * 10}")
            else:
                logging.debug("Multithreads map:")
                logging.debug(f"{'-' * 10} Source grid {'-' * 10}")
                for key in ['src_global_x', 'src_global_y', 'src_global_x_size', 'src_global_y_size', ]:
                    logging.debug(f"{key} = {self.parallel_map[self.rank][key]}")
                logging.debug(f"{'-' * 10} Target grid {'-' * 10}")
                for key in ['dst_global_x', 'dst_global_y', 'dst_local_x_size', 'dst_local_y_size']:
                    logging.debug(f"{key} = {self.parallel_map[self.rank][key]}")

            for key in self.parallel_map[self.rank]:
                if MPI_FOUND and key != "threads":
                    logging.debug(f"{key} = {self.parallel_map[self.rank][key]}")

            # if len(self.parallel_map[self.rank]["threads"]) != 1:
            for thread in range(len(self.parallel_map[self.rank]["threads"])):
                logging.debug(f"   {'-' * 10} Rank {self.rank} - Thread n° {thread} {'-' * 10}")
                for thread_key in self.parallel_map[self.rank]["threads"][thread]:
                    logging.debug(
                        f"    {thread_key} = {self.parallel_map[self.rank]['threads'][thread][thread_key]}")

            if self.rank == 0:
                logging.debug("-" * 20)

            if self.horizontal_resampling:
                self.compute_weight()

        # try to fill metadata
        self.read_metadata()

    def check_bbox_validity(self, candidate):
        """
        Check the validity of the bounding box.

        Parameters
        ----------
        candidate : list
            Bounding box coordinates [xmin, xmax, ymin, ymax].

        Returns
        -------
        bool
            True if the bounding box is valid, False otherwise.

        Examples
        --------
        >>> coverage = Coverage(myReader, bbox=[0, 10, 0, 10], resolution_x=1.0, resolution_y=1.0)
        """
        # TODO improve error messages

        Ymin = candidate[2]
        Ymax = candidate[3]
        Xmin = candidate[0]
        Xmax = candidate[1]

        if Xmax <= Xmin:
            return False
        if Ymax <= Ymin:
            return False

        if Ymin < np.min(self.source_global_axis_y):
            return False
        if Ymax > np.max(self.source_global_axis_y):
            return False
        if Xmin < np.min(self.source_global_axis_x):
            return False
        if Xmax > np.max(self.source_global_axis_x):
            return False

        return True

    def check_point_is_inside(self, target_lon, target_lat, lon, lat, tolerance=5):
        """
        Check if a point is inside the coverage.

        Parameters
        ----------
        target_lon : float
            Longitude of the target point.
        target_lat : float
            Latitude of the target point.
        lon : array
            Longitude array.
        lat : array
            Latitude array.
        tolerance : int, optional
            Tolerance for rounding.

        Returns
        -------
        bool
            True if the point is inside the coverage, False otherwise.

        Examples
        --------
        >>> coverage.check_point_is_inside(5.0, 45.0, lon_array, lat_array)
        True
        """

        if np.round(target_lat, decimals=tolerance) < np.min(np.round(lat, decimals=tolerance)):
            return False
        if np.round(target_lat, decimals=tolerance) > np.max(np.round(lat, decimals=tolerance)):
            return False
        if np.round(target_lon, decimals=tolerance) < np.min(np.round(lon, decimals=tolerance)):
            return False
        if np.round(target_lon, decimals=tolerance) > np.max(np.round(lon, decimals=tolerance)):
            return False

        return True

    def __init_parallel_map(self):
        """
        Create the parallel map for parallel processing.
        The parallel map is a dictionary that contains array slice of the source and destination grids for each MPI rank and sub threads.

        Example of 4 x 4 grid with MPI only

        Source global grid with the MPI splitting
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+

        MPI Rank n°0
        src_global_x = slice(0, 4, None) ; array slice for the x axis in the global source grid
        src_global_y = slice(0, 2, None) ; array slice for the y axis in the global source grid
        src_global_x_size = 4 ; size of the array slice for the x axis in the global source grid
        src_global_y_size = 2 ; size of the array slice for the y axis in the global source grid
        src_global_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the global source grid
        src_global_y_overlap = slice(0, 4, None)  ; array slice with overlapping for the y axis in the global source grid
        src_global_x_size_overlap = 4  ; size of the array slice with overlapping for the x axis in the global source grid
        src_global_y_size_overlap = 4 ; size of the array slice with overlapping for the y axis in the global source grid

        Source local grid of MPI rank n°0
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        overlap with rank n°1
        +----+----+----+----+

        src_local_x = slice(0, 4, None) ; array slice for the x axis in the local source grid
        src_local_y = slice(0, 2, None) ; array slice for the y axis in the local source grid
        src_local_x_size = 4 ; size of the array slice for the x axis in the local source grid
        src_local_y_size = 2 ; size of the array slice for the y axis in the local source grid
        src_local_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the local source grid
        src_local_y_overlap = slice(0, 4, None) ; array slice with overlapping for the y axis in the local source grid
        src_local_x_size_overlap = 4  ; size of the array slice with overlapping for the x axis in the local source grid
        src_local_y_size_overlap = 4  ; size of the array slice with overlapping for the y axis in the local source grid

        MPI Rank n°1
        src_global_x = slice(0, 4, None)     ; array slice for the x axis in the global source grid
        src_global_y = slice(2, 4, None)     ; array slice for the y axis in the global source grid
        src_global_x_size = 4                ; size of the array slice for the x axis in the global source grid
        src_global_y_size = 2                ; size of the array slice for the y axis in the global source grid
        src_global_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the global source grid
        src_global_y_overlap = slice(0, 4, None) ; array slice with overlapping for the y axis in the global source grid
        src_global_x_size_overlap = 4    ; size of the array slice with overlapping for the x axis in the global source grid
        src_global_y_size_overlap = 4    ; size of the array slice with overlapping for the y axis in the global source grid

        Source local grid of MPI rank n°1
        +----+----+----+----+
        overlap with rank n°0
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+

        src_local_x = slice(0, 4, None) ; array slice for the x axis in the local source grid
        src_local_y = slice(2, 4, None) ; array slice for the y axis in the local source grid
        src_local_x_size = 4 ; size of the array slice for the x axis in the local source grid
        src_local_y_size = 2 ; size of the array slice for the y axis in the local source grid
        src_local_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the local source grid
        src_local_y_overlap = slice(0, 4, None) ; array slice with overlapping for the y axis in the local source grid
        src_local_x_size_overlap = 4  ; size of the array slice with overlapping for the x axis in the local source grid
        src_local_y_size_overlap = 4  ; size of the array slice with overlapping for the y axis in the local source grid

        Destination global grid
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+

        MPI Rank n°0
        dst_global_x = slice(0, 4, None) ; array slice for the x axis in the global destination grid
        dst_global_y = slice(0, 2, None) ; array slice for the y axis in the global destination grid
        dst_global_x_size = 4 ; size of the array slice for the x axis in the global destination grid
        dst_global_y_size = 2 ; size of the array slice for the y axis in the global destination grid
        dst_global_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the global destination grid
        dst_global_y_overlap = slice(0, 4, None)  ; array slice with overlapping for the y axis in the global destination grid
        dst_global_x_size_overlap = 4  ; size of the array slice with overlapping for the x axis in the global destination grid
        dst_global_y_size_overlap = 4 ; size of the array slice with overlapping for the y axis in the global destination grid

        Destination local grid of MPI rank n°0
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        | 0  | 0  | 0  | 0  |
        +----+----+----+----+
        overlap with rank n°1
        +----+----+----+----+

        dst_local_x = slice(0, 4, None) ; array slice for the x axis in the local destination grid
        dst_local_y = slice(0, 2, None) ; array slice for the y axis in the local destination grid
        dst_local_x_size = 4 ; size of the array slice for the x axis in the local destination grid
        dst_local_y_size = 2 ; size of the array slice for the y axis in the local destination grid
        dst_local_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the local destination grid
        dst_local_y_overlap = slice(0, 4, None) ; array slice with overlapping for the y axis in the local destination grid
        dst_local_x_size_overlap = 4  ; size of the array slice with overlapping for the x axis in the local destination grid
        dst_local_y_size_overlap = 4  ; size of the array slice with overlapping for the y axis in the local destination grid

        MPI Rank n°1
        dst_global_x = slice(0, 4, None)     ; array slice for the x axis in the global destination grid
        dst_global_y = slice(2, 4, None)     ; array slice for the y axis in the global destination grid
        dst_global_x_size = 4                ; size of the array slice for the x axis in the global destination grid
        dst_global_y_size = 2                ; size of the array slice for the y axis in the global destination grid
        dst_global_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the global destination grid
        dst_global_y_overlap = slice(0, 4, None) ; array slice with overlapping for the y axis in the global destination grid
        dst_global_x_size_overlap = 4    ; size of the array slice with overlapping for the x axis in the global destination grid
        dst_global_y_size_overlap = 4    ; size of the array slice with overlapping for the y axis in the global destination grid

        Destination local grid of MPI rank n°1
        +----+----+----+----+
        overlap with rank n°0
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+
        | 1  | 1  | 1  | 1  |
        +----+----+----+----+

        dst_local_x = slice(0, 4, None) ; array slice for the x axis in the local destination grid
        dst_local_y = slice(2, 4, None) ; array slice for the y axis in the local destination grid
        dst_local_x_size = 4 ; size of the array slice for the x axis in the local destination grid
        dst_local_y_size = 2 ; size of the array slice for the y axis in the local destination grid
        dst_local_x_overlap = slice(0, 4, None) ; array slice with overlapping for the x axis in the local destination grid
        dst_local_y_overlap = slice(0, 4, None) ; array slice with overlapping for the y axis in the local destination grid
        dst_local_x_size_overlap = 4  ; size of the array slice with overlapping for the x axis in the local destination grid
        dst_local_y_size_overlap = 4  ; size of the array slice with overlapping for the y axis in the local destination grid

        Examples
        --------
        >>> coverage.__init_parellel_map()
        """
        target_mpi_sample = (self.target_global_y_size, self.target_global_x_size)

        # Split the axes with the MPI size
        target_mpi_slices = shape_split(target_mpi_sample, self.size, axis=[0, 0])
        # If we can divide the grid dimensions with the number of MPI rank,
        # we raise an Error
        if self.size != len(target_mpi_slices.flatten()):
            logging.error(
                f"Unable to divide the grid {target_mpi_sample} with the MPI size ({self.size}). Try with different MPI size.")
            raise ValueError(
                f"Unable to divide the grid {target_mpi_sample} with the MPI size ({self.size}). Try with different MPI size.")

        mpi_slice_index = 0
        for slyce in target_mpi_slices.flatten():
            mpi_slice = tuple(slyce)
            self.parallel_map[mpi_slice_index] = self.compute_slice_coordinates(
                mpi_slice,
                self.source_global_x_size,
                self.source_global_y_size,
                self.target_global_x_size,
                self.target_global_y_size)

            # Split the axes with the number of threads
            target_threads_sample = (self.parallel_map[mpi_slice_index]["dst_local_y_size"],
                                     self.parallel_map[mpi_slice_index]["dst_local_x_size"])
            target_threads_slices = shape_split(target_threads_sample, self.threads_number, axis=[0, 0])

            # If we can divide the grid dimensions with the number of threads,
            # we set the threads number with the number of slice
            self.threads_number = len(target_threads_slices.flatten())
            logging.debug(f"Threads number is {self.threads_number}")

            self.parallel_map[mpi_slice_index]['threads'] = np.empty([self.threads_number], dtype=object)

            slice_thread_index = 0
            for thread_slyce in target_threads_slices.flatten():
                thread_slice = tuple(thread_slyce)
                self.parallel_map[mpi_slice_index]["threads"][slice_thread_index] = self.compute_slice_coordinates(
                    thread_slice,
                    self.parallel_map[mpi_slice_index]["src_local_x_size"],
                    self.parallel_map[mpi_slice_index]["src_local_y_size"],
                    self.parallel_map[mpi_slice_index]["dst_local_x_size"],
                    self.parallel_map[mpi_slice_index]["dst_local_y_size"],
                    self.parallel_map[mpi_slice_index]
                )
                slice_thread_index = slice_thread_index + 1

            mpi_slice_index = mpi_slice_index + 1

    def compute_slice_coordinates(self,
                                  slice: tuple,
                                  source_global_x_size: int,
                                  source_global_y_size: int,
                                  target_global_x_size: int,
                                  target_global_y_size: int,
                                  parent_slice: tuple = None):
        """
        Compute slice coordinates in the source grid and the destination grid with overlap.
        Args:
            slice (tuple): Current slice to compute
            source_global_x_size (int): Size of the global source x axis
            source_global_y_size (int): Size of the global source y axis
            target_global_x_size (int): Size of the global target x axis
            target_global_y_size (int): Size of the global target y axis
            parent_slice (tuple, optional): Slice of the parent slice
        """

        map = {}

        ### Destination grid ###
        if parent_slice is not None:
            map["dst_global_x"] = np.s_[
                parent_slice["dst_global_x"].start + slice[1].start:
                parent_slice["dst_global_x"].start + slice[1].stop
            ]

            map["dst_global_y"] = np.s_[
                parent_slice["dst_global_y"].start + slice[0].start:
                parent_slice["dst_global_y"].start + slice[0].stop
            ]
        else:
            map["dst_global_x"] = slice[1]
            map["dst_global_y"] = slice[0]

        # Compute overlap
        if parent_slice is not None:

            # Position of the thread core inside the parent overlap grid
            thread_x_start = parent_slice["dst_local_x"].start + slice[1].start
            thread_x_stop = parent_slice["dst_local_x"].start + slice[1].stop

            thread_y_start = parent_slice["dst_local_y"].start + slice[0].start
            thread_y_stop = parent_slice["dst_local_y"].start + slice[0].stop

            # X
            dst_parent_x_min_overlap = max(
                0,
                thread_x_start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_parent_x_max_overlap = min(
                parent_slice["dst_global_x_size_overlap"],
                thread_x_stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_global_x_min_overlap = max(
                parent_slice["dst_global_x_overlap"].start,
                map["dst_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_global_x_max_overlap = min(
                parent_slice["dst_global_x_overlap"].stop,
                map["dst_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            # Y
            dst_parent_y_min_overlap = max(
                0,
                thread_y_start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_parent_y_max_overlap = min(
                parent_slice["dst_global_y_size_overlap"],
                thread_y_stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_global_y_min_overlap = max(
                parent_slice["dst_global_y_overlap"].start,
                map["dst_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_global_y_max_overlap = min(
                parent_slice["dst_global_y_overlap"].stop,
                map["dst_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

        else:

            # X
            dst_global_x_min_overlap = max(
                0,
                map["dst_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_global_x_max_overlap = min(
                target_global_x_size,
                map["dst_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            # Y
            dst_global_y_min_overlap = max(
                0,
                map["dst_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            dst_global_y_max_overlap = min(
                target_global_y_size,
                map["dst_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

        map["dst_global_x_overlap"] = np.s_[
            dst_global_x_min_overlap:dst_global_x_max_overlap
        ]

        map["dst_global_y_overlap"] = np.s_[
            dst_global_y_min_overlap:dst_global_y_max_overlap
        ]

        map["dst_global_x_size_overlap"] = (
                map["dst_global_x_overlap"].stop -
                map["dst_global_x_overlap"].start
        )

        map["dst_global_y_size_overlap"] = (
                map["dst_global_y_overlap"].stop -
                map["dst_global_y_overlap"].start
        )

        if parent_slice is not None:
            map["dst_parent_x"] = np.s_[thread_x_start:thread_x_stop]
            map["dst_parent_y"] = np.s_[thread_y_start:thread_y_stop]

            map["dst_parent_x_overlap"] = np.s_[
                dst_parent_x_min_overlap:dst_parent_x_max_overlap
            ]

            map["dst_parent_y_overlap"] = np.s_[
                dst_parent_y_min_overlap:dst_parent_y_max_overlap
            ]

            map["dst_parent_x_size_overlap"] = (
                    map["dst_parent_x_overlap"].stop -
                    map["dst_parent_x_overlap"].start
            )

            map["dst_parent_y_size_overlap"] = (
                    map["dst_parent_y_overlap"].stop -
                    map["dst_parent_y_overlap"].start
            )

        # Compute destination local grid
        dst_local_x_size = (
                map["dst_global_x"].stop -
                map["dst_global_x"].start
        )

        dst_local_y_size = (
                map["dst_global_y"].stop -
                map["dst_global_y"].start
        )

        if map["dst_global_x"].start <= Coverage.HORIZONTAL_OVERLAPING_SIZE:
            dst_local_x_min = map["dst_global_x"].start
        else:
            dst_local_x_min = (
                    map["dst_global_x"].start -
                    map["dst_global_x_overlap"].start
            )

        if map["dst_global_y"].start <= Coverage.HORIZONTAL_OVERLAPING_SIZE:
            dst_local_y_min = map["dst_global_y"].start
        else:
            dst_local_y_min = (
                    map["dst_global_y"].start -
                    map["dst_global_y_overlap"].start
            )

        dst_local_x_max = dst_local_x_min + dst_local_x_size
        dst_local_y_max = dst_local_y_min + dst_local_y_size

        map["dst_local_x"] = np.s_[
            dst_local_x_min:dst_local_x_max
        ]

        map["dst_local_y"] = np.s_[
            dst_local_y_min:dst_local_y_max
        ]

        map["dst_local_x_size"] = dst_local_x_size
        map["dst_local_y_size"] = dst_local_y_size

        dst_local_x_min_overlap = max(
            0,
            map["dst_local_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
        )

        dst_local_x_max_overlap = min(
            target_global_x_size,
            map["dst_local_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
        )

        dst_local_y_min_overlap = max(
            0,
            map["dst_local_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
        )

        dst_local_y_max_overlap = min(
            target_global_y_size,
            map["dst_local_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
        )

        map["dst_local_x_overlap"] = np.s_[
            dst_local_x_min_overlap:dst_local_x_max_overlap
        ]

        map["dst_local_y_overlap"] = np.s_[
            dst_local_y_min_overlap:dst_local_y_max_overlap
        ]

        map["dst_local_x_size_overlap"] = (
                map["dst_local_x_overlap"].stop -
                map["dst_local_x_overlap"].start
        )

        map["dst_local_y_size_overlap"] = (
                map["dst_local_y_overlap"].stop -
                map["dst_local_y_overlap"].start
        )

        # TODO case of target irregular grid
        target_global_axis_x = self.target_global_axis_x[
            map["dst_global_x"]
        ]

        target_global_axis_y = self.target_global_axis_y[
            map["dst_global_y"]
        ]

        ### Source grid ###
        # Find xmin, xmax, ymin, ymax coordinates
        # of the destination grid in the source grid
        if self.is_regular_grid(type="source"):

            if np.min(target_global_axis_x) == np.max(target_global_axis_x):

                idx = np.where(
                    self.source_global_axis_x >= np.min(target_global_axis_x)
                )

            else:

                idx = np.where(
                    (self.source_global_axis_x >= np.min(target_global_axis_x)) &
                    (self.source_global_axis_x <= np.max(target_global_axis_x))
                )

            if idx[0].size == 0:
                raise ValueError(
                    "Unable to find target xmin / xmax in the source grid"
                )

            xmin = np.min(idx[0])
            xmax = np.max(idx[0]) + 1

            if np.min(target_global_axis_y) == np.max(target_global_axis_y):

                idx = np.where(
                    self.source_global_axis_y >= np.min(target_global_axis_y)
                )

            else:

                idx = np.where(
                    (self.source_global_axis_y >= np.min(target_global_axis_y)) &
                    (self.source_global_axis_y < np.max(target_global_axis_y))
                )

            if idx[0].size == 0:
                raise ValueError(
                    "Unable to find target ymin / ymax in the source grid"
                )

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1

        else:

            idx = np.where(
                (self.source_global_axis_x >= np.min(target_global_axis_x)) &
                (self.source_global_axis_x <= np.max(target_global_axis_x)) &
                (self.source_global_axis_y >= np.min(target_global_axis_y)) &
                (self.source_global_axis_y <= np.max(target_global_axis_y))
            )

            if np.shape(idx)[1] == 0:

                logging.debug(
                    "Unable to find target xmin / xmax / ymin / ymax in the source grid"
                )

                ymin = 0
                ymax = self.source_global_y_size
                xmin = 0
                xmax = self.source_global_x_size

            else:

                ymin = np.min(idx[0])
                ymax = np.max(idx[0]) + 1

                xmin = np.min(idx[1])
                xmax = np.max(idx[1]) + 1

        # Compute source global slice.
        # xmin/xmax/ymin/ymax are already global indices because
        # they were computed from source_global_axis_x/y.
        map["src_global_x"] = np.s_[
            int(xmin):int(xmax)
        ]

        map["src_global_y"] = np.s_[
            int(ymin):int(ymax)
        ]

        map["src_global_x_size"] = xmax - xmin
        map["src_global_y_size"] = ymax - ymin

        # --------------------------------------------------
        # Source global X overlap
        # --------------------------------------------------
        if parent_slice is not None:

            parent_src_x_overlap = parent_slice["src_global_x_overlap"]

            src_global_x_min_overlap = max(
                parent_src_x_overlap.start,
                map["src_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            src_global_x_max_overlap = min(
                parent_src_x_overlap.stop,
                map["src_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

        else:

            src_global_x_min_overlap = max(
                0,
                map["src_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            src_global_x_max_overlap = min(
                source_global_x_size,
                map["src_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

        map["src_global_x_overlap"] = np.s_[
            src_global_x_min_overlap:src_global_x_max_overlap
        ]

        # --------------------------------------------------
        # Source global Y overlap
        # --------------------------------------------------
        if parent_slice is not None:

            parent_src_y_overlap = parent_slice["src_global_y_overlap"]

            src_global_y_min_overlap = max(
                parent_src_y_overlap.start,
                map["src_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            src_global_y_max_overlap = min(
                parent_src_y_overlap.stop,
                map["src_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

        else:

            src_global_y_min_overlap = max(
                0,
                map["src_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

            src_global_y_max_overlap = min(
                source_global_y_size,
                map["src_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE
            )

        map["src_global_y_overlap"] = np.s_[
            src_global_y_min_overlap:src_global_y_max_overlap
        ]

        # Source global overlap size
        map["src_global_x_size_overlap"] = (
                map["src_global_x_overlap"].stop -
                map["src_global_x_overlap"].start
        )

        map["src_global_y_size_overlap"] = (
                map["src_global_y_overlap"].stop -
                map["src_global_y_overlap"].start
        )

        # Recompute source local slice
        src_local_x_size = map["src_global_x_size"]
        src_local_y_size = map["src_global_y_size"]

        if map["src_global_x"].start <= Coverage.HORIZONTAL_OVERLAPING_SIZE:
            src_local_x_min = map["src_global_x"].start
        else:
            src_local_x_min = (
                    map["src_global_x"].start -
                    map["src_global_x_overlap"].start
            )

        src_local_x_max = src_local_x_min + src_local_x_size

        if map["src_global_y"].start <= Coverage.HORIZONTAL_OVERLAPING_SIZE:
            src_local_y_min = map["src_global_y"].start
        else:
            src_local_y_min = (
                    map["src_global_y"].start -
                    map["src_global_y_overlap"].start
            )

        src_local_y_max = src_local_y_min + src_local_y_size

        map["src_local_x"] = np.s_[
            int(src_local_x_min):int(src_local_x_max)
        ]

        map["src_local_y"] = np.s_[
            int(src_local_y_min):int(src_local_y_max)
        ]

        map["src_local_x_size"] = src_local_x_size
        map["src_local_y_size"] = src_local_y_size

        # Overlap
        src_local_x_size_overlap = map["src_global_x_size_overlap"]
        src_local_y_size_overlap = map["src_global_y_size_overlap"]

        map["src_local_x_overlap"] = np.s_[
            0:src_local_x_size_overlap
        ]

        map["src_local_y_overlap"] = np.s_[
            0:src_local_y_size_overlap
        ]

        map["src_local_x_size_overlap"] = src_local_x_size_overlap
        map["src_local_y_size_overlap"] = src_local_y_size_overlap

        return map

    def compute_weight(self):
        """
        Compute horizontal interpolation weight (Delaunay triangulation)
        """
        if self.horizontal_resampling:
            self.source_global_tri = np.empty([self.size], dtype=object)

            if self.rank == 0:
                logging.info(
                    '[horizontal_interpolation] Compute weights...')

            self.source_global_tri[self.rank] = interp_2d_weights(
                self.read_axis_x(type="source", with_overlap=True, ),
                self.read_axis_y(type="source", with_overlap=True),
                self.parallel_map[self.rank]["threads"],
                Coverage.HORIZONTAL_INTERPOLATOR,
            )

    def read_metadata(self):
        """
        Read the metadata from the file if the reader contains a read_metadata() function.
        Examples
        --------
        >>> coverage.read_metadata()
        """
        if "read_metadata" in dir(self.reader):
            m = self.reader.read_metadata()

    def get_x_size(self, type="target", with_overlap=False):
        """
        Get the x size.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        int
            The x size.
        Examples
        --------
        >>> coverage.get_x_size()
        120
        >>> coverage.get_x_size(type="source", with_overlap=True)
        124
        """
        if type == "target_global":
            return self.target_global_x_size
        elif type == "source_global":
            return self.source_global_x_size
        elif type == "source" and with_overlap is True:
            return self.parallel_map[self.rank]["src_local_x_size_overlap"]
        elif type == "source" and with_overlap is False:
            return self.parallel_map[self.rank]["src_local_x_size"]
        elif type == "target" and with_overlap is True:
            return self.parallel_map[self.rank]["dst_local_x_size_overlap"]
        else:
            return self.parallel_map[self.rank]["dst_local_x_size"]

    def get_y_size(self, type="target", with_overlap=False):
        """
        Get the y size.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        int
            The y size.

        Examples
        --------
        >>> coverage.get_y_size()
        80
        >>> coverage.get_y_size(type="source", with_overlap=True)
        84
        """
        if type == "target_global":
            return self.target_global_y_size
        elif type == "source_global":
            return self.source_global_y_size
        elif type == "source" and with_overlap is True:
            return self.parallel_map[self.rank]["src_local_y_size_overlap"]
        elif type == "source" and with_overlap is False:
            return self.parallel_map[self.rank]["src_local_y_size"]
        elif type == "target" and with_overlap is True:
            return self.parallel_map[self.rank]["dst_local_y_size_overlap"]
        else:
            return self.parallel_map[self.rank]["dst_local_y_size"]

    def is_regular_grid(self, type="target"):
        """
        Check if the grid is regular.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source").

        Returns
        -------
        bool
            True if the grid is regular, False otherwise.
        Examples
        --------
        >>> coverage.is_regular_grid()
        True
        >>> coverage.is_regular_grid(type="source")
        False
        """
        if type == "target":
            return self.target_regular_grid
        else:
            return self.source_regular_grid

    # Axis
    def read_axis_x(self, type="target", with_overlap=False):
        """
        Return the values of the x axis.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        array
            A one- or two-dimensional array of x axis values (often longitude): [x] or [y, x], depending on the grid type.

        Examples
        --------
        >>> x_axis = coverage.read_axis_x()
        >>> x_axis = coverage.read_axis_x(type="source", with_overlap=True)
        """

        if type == "target_global":
            return self.target_global_axis_x

        elif type == "source_global":
            return self.source_global_axis_x

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_x(self.parallel_map[self.rank]["src_global_x_overlap"].start,
                                           self.parallel_map[self.rank]["src_global_x_overlap"].stop,
                                           self.parallel_map[self.rank]["src_global_y_overlap"].start,
                                           self.parallel_map[self.rank]["src_global_y_overlap"].stop)
        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_x(self.parallel_map[self.rank]["src_global_x"].start,
                                           self.parallel_map[self.rank]["src_global_x"].stop,
                                           self.parallel_map[self.rank]["src_global_y"].start,
                                           self.parallel_map[self.rank]["src_global_y"].stop)
        elif type == "target" and with_overlap is True:
            if self.is_regular_grid():
                return self.target_global_axis_x[self.parallel_map[self.rank]["dst_global_x_overlap"]]
            else:
                return self.target_global_axis_x[self.parallel_map[self.rank]["dst_global_y_overlap"],
                self.parallel_map[self.rank]["dst_global_x_overlap"]]
        elif type == "target" and with_overlap is False:
            if self.is_regular_grid():
                return self.target_global_axis_x[self.parallel_map[self.rank]["dst_global_x"]]
            else:
                return self.target_global_axis_x[
                    self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]

    def read_axis_y(self, type="target", with_overlap=False):
        """
        Return the values (often latitude) of the y axis.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        array
            A one- or two-dimensional array of y axis values (often latitude): [x] or [y, x], depending on the grid type.

        Examples
        --------
        >>> y_axis = coverage.read_axis_y()
        >>> y_axis = coverage.read_axis_y(type="source", with_overlap=True)
        """
        if type == "target_global":
            return self.target_global_axis_y

        elif type == "source_global":
            return self.source_global_axis_y

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_y(self.parallel_map[self.rank]["src_global_x_overlap"].start,
                                           self.parallel_map[self.rank]["src_global_x_overlap"].stop,
                                           self.parallel_map[self.rank]["src_global_y_overlap"].start,
                                           self.parallel_map[self.rank]["src_global_y_overlap"].stop)
        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_y(self.parallel_map[self.rank]["src_global_x"].start,
                                           self.parallel_map[self.rank]["src_global_x"].stop,
                                           self.parallel_map[self.rank]["src_global_y"].start,
                                           self.parallel_map[self.rank]["src_global_y"].stop)
        elif type == "target" and with_overlap is True:
            if self.is_regular_grid():
                return self.target_global_axis_y[self.parallel_map[self.rank]["dst_global_y_overlap"]]
            else:
                return self.target_global_axis_y[self.parallel_map[self.rank]["dst_global_y_overlap"],
                self.parallel_map[self.rank]["dst_global_x_overlap"]]
        elif type == "target" and with_overlap is False:
            if self.is_regular_grid():
                return self.target_global_axis_y[self.parallel_map[self.rank]["dst_global_y"]]
            else:
                return self.target_global_axis_y[
                    self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]

    def find_point_index(self, target_lon, target_lat, decimal_tolerance=5, method="classic", only_mask_value=True,
                         type="source"):
        """
                Find the index of the closest point to the given point.

                Parameters
                ----------
                target_lon : float
                    Longitude of the target point.
                target_lat : float
                    Latitude of the target point.
                decimal_tolerance : int, optional
                    Tolerance for rounding coordinates.
                method : str, optional
                    Method for calculation ("classic", "quick").
                only_mask_value : bool, optional
                    Whether to consider only mask values.
                type : str, optional
                    Type of the grid ("source", "source_global").

                Returns
                -------
                list
                    A list containing:
                    [0] : the x index of the closest point
                    [1] : the y index of the closest point
                    [2] : the longitude of the closest point
                    [3] : the latitude of the closest point
                    [4] : the distance of the closest point in kilometers.

                Examples
                --------
                >>> idx = coverage.find_point_index(5.0, 45.0)
                >>> idx = coverage.find_point_index(5.0, 45.0, method="quick", type="source_global")
                """
        lon = self.read_axis_x(type="source", with_overlap=False)
        lat = self.read_axis_y(type="source", with_overlap=False)

        if self.check_point_is_inside(target_lon, target_lat, lon, lat, tolerance=decimal_tolerance):
            try:
                mask = self.read_variable_2D_sea_binary_mask(type="source", with_overlap=False)
            except NotImplementedError:
                logging.warning("No 2D sea binary mask found")
                # mask = np.ones([self.source_global_y_size, self.source_global_x_size])
                only_mask_value = False

            dist = np.zeros([self.get_y_size(type="source", with_overlap=False),
                             self.get_x_size(type="source", with_overlap=False)])
            dist[:] = 10000000

            if method == "classic":
                for x in range(0, self.get_x_size(type="source", with_overlap=False)):
                    for y in range(0, self.get_y_size(type="source", with_overlap=False)):

                        if only_mask_value:
                            if (mask[y, x] == 1):  # =Terre
                                if self.reader.is_regular_grid():
                                    dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[x], lat[y])
                                else:
                                    dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[y, x], lat[y, x])
                        else:
                            if self.reader.is_regular_grid():
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[x], lat[y])
                            else:
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[y, x], lat[y, x])

                nearest_y_index, nearest_x_index = np.where(dist == np.min(dist))

                if len(nearest_y_index) == 0 or len(nearest_x_index) == 0:
                    logging.error("no neaarest point found")
                    raise RuntimeError("No nearest point found")

                nearest_x_index = nearest_x_index[0]
                nearest_y_index = nearest_y_index[0]
                min_dist = dist[nearest_y_index, nearest_x_index]

                if self.is_regular_grid(type="source"):
                    nearest_lon = lon[nearest_x_index]
                    nearest_lat = lat[nearest_y_index]
                else:
                    nearest_lon = lon[nearest_y_index, nearest_x_index]
                    nearest_lat = lat[nearest_y_index, nearest_x_index]

                if type == "source":
                    return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]
                elif type == "source_global":
                    return [self.parallel_map[self.rank]["src_global_x"].start + nearest_x_index,
                            self.parallel_map[self.rank]["src_global_y"].start + nearest_y_index, nearest_lon,
                            nearest_lat,
                            min_dist]
                else:
                    raise ValueError("Type doesn't match [source, source_global]")

            elif method == "quick":

                if self.is_regular_grid():

                    # Longitude : on cherche l'index le plus proche
                    array = np.asarray(lon)
                    nearest_x_index = (np.abs(array - target_lon)).argmin()

                    # Latitude : on cherche l'index le plus proche
                    array = np.asarray(lat)
                    nearest_y_index = (np.abs(array - target_lat)).argmin()

                    nearest_lon = lon[nearest_x_index]
                    nearest_lat = lat[nearest_y_index]

                    min_dist = distance_on_unit_sphere(target_lon, target_lat, nearest_lon, nearest_lat)

                    if type == "source":
                        return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]
                    elif type == "source_global":
                        return [self.parallel_map[self.rank]["src_global_x"].start + nearest_x_index,
                                self.parallel_map[self.rank]["src_global_y"].start + nearest_y_index, nearest_lon,
                                nearest_lat, min_dist]
                    else:
                        raise ValueError("Type doesn't match [source, source_global]")

                else:
                    raise NotImplementedError("Method " + str(method) + " is not implemented for regular grid.")

            else:
                raise RuntimeError("Method " + str(method) + " is not implemented yet.")
        else:
            logging.warning("Point is outside the rank n°" + str(self.rank))
            raise NotFoundInRankError(self.rank, "Point is outside the rank")

    def resample_2d_variable(self, values, is_binary_data=False):
        return resample_2d_to_grid(
            self.source_global_tri[self.rank],
            self.read_axis_x(type="target", with_overlap=True),
            self.read_axis_y(type="target", with_overlap=True),
            values,
            Coverage.HORIZONTAL_INTERPOLATION_METHOD,
            self.parallel_map[self.rank]["threads"],
            Coverage.HORIZONTAL_INTERPOLATOR,
            is_binary_data
        )[self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    def __read_variable(self, function_name, is_binary_data=False):

        fn = getattr(self.reader, function_name)

        data = fn(
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        is_vector = True if len(np.shape(data)) == 3 and np.shape(data)[0] == 2 else False

        if self.horizontal_resampling:
            if is_vector:
                return [
                    self.resample_2d_variable(data[0], is_binary_data=is_binary_data),
                    self.resample_2d_variable(data[1], is_binary_data=is_binary_data)
                ]
            else:
                return self.resample_2d_variable(data, is_binary_data=is_binary_data)
        else:
            if is_vector:
                return [
                    data[0][self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]],
                    data[1][self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]
                ]
            else:
                return data[self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]

    # Variables
    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self):
        """
        Read the bathymetry over the entire coverage.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> bathymetry = coverage.read_variable_bathymetry()
        """
        return self.__read_variable(inspect.stack()[0][3])

    def read_variable_topography(self, ):
        """
        Read the topography over the entire coverage.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> topography = coverage.read_variable_topography()
        """
        return self.__read_variable(inspect.stack()[0][3])

    def read_variable_mesh_size(self):
        """
        Read the mesh size over the entire coverage.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> mesh_size = coverage.read_variable_mesh_size()
        """
        return self.__read_variable(inspect.stack()[0][3])

    def read_variable_x_mesh_size(self):
        """
        Read the mesh size on the x axis.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> x_mesh_size = coverage.read_variable_x_mesh_size()
        """
        return self.__read_variable(inspect.stack()[0][3])

    def read_variable_y_mesh_size(self):
        """
        Read the mesh size on the y axis.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> y_mesh_size = coverage.read_variable_y_mesh_size()
        """
        return self.__read_variable(inspect.stack()[0][3])

    def read_variable_2D_sea_binary_mask(self):
        """
        Read the land/sea mask over the entire coverage.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        array
            A two-dimensional array [y, x].
            0 = Land
            1 = Sea
        Examples
        --------
        >>> mask = coverage.read_variable_2D_sea_binary_mask()
        """
        return self.__read_variable(inspect.stack()[0][3], is_binary_data=True)

    def read_variable_Ha(self):
        """
        Read tide amplitude

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> Ha = coverage.read_variable_Ha()
        """
        return self.__read_variable(inspect.stack()[0][3])
