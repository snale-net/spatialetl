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
import pytest

from spatialetl.coverage.coverage import Coverage
from spatialetl.coverage.io.memory_reader import MemoryReader
from spatialetl.utils.logger import logging

# Parameters
MATRIX_SIZES = [5, 10, 50]
RESOLUTIONS = [1.0, 0.5, 0.05]
THREADS = [1, 2, 3, 4, 6, 8, 10, 12]
OVERLAPS = [2, 4, 7, 10]

INTERPOLATORS = [
    "cgal",
    "scipy",
]

@pytest.fixture
def test_data(request):
    size, resolution, affine = request.param

    # Source grid
    axis = np.arange(size, dtype=float)
    source_x_axis, source_y_axis = np.meshgrid(axis, axis)

    if affine:
        data = 2 * source_x_axis + 3 * source_y_axis
        data_function = lambda x, y: 2 * x + 3 * y
    else:
        data_function = lambda x, y: (
                np.sin(x / 2.0)
                + np.cos(y / 3.0)
                + 0.2 * np.sin((x + y) / 2.0)
                + 0.01 * x * y
        )
        data = data_function(source_x_axis, source_y_axis)

    # Target grid
    target_x_axis = np.arange(0, size - 1, resolution)
    target_y_axis = np.arange(0, size - 1, resolution)

    target_x, target_y = np.meshgrid(
        target_x_axis,
        target_y_axis
    )

    target_data = data_function(target_x, target_y)

    return {
        "source_x": source_x_axis,
        "source_y": source_y_axis,
        "data": data,
        "target_x": target_x_axis,
        "target_y": target_y_axis,
        "target_data": target_data,
        "resolution": resolution,
    }


@pytest.mark.parametrize(
    "interpolator",
    INTERPOLATORS
)
@pytest.mark.parametrize(
    "test_data",
    [
        (size, resolution, True)  # Test with affine data
        for size in MATRIX_SIZES
        for resolution in RESOLUTIONS
    ],
    indirect=True
)
def test_sequential_interpolation(test_data, interpolator, monkeypatch, caplog):
    caplog.set_level(logging.INFO)

    logging.info(
        f"Testing size={test_data['source_x'].shape}, "
        f"resolution={test_data['resolution']}, "
        f"interpolator={interpolator}, "
    )

    monkeypatch.setattr(
        Coverage,
        "HORIZONTAL_INTERPOLATOR",
        interpolator
    )

    reader = MemoryReader(
        x_axis=test_data["source_x"],
        y_axis=test_data["source_y"],
        bathymetry=test_data["data"]
    )

    coverage = Coverage(
        reader=reader,
        resolution_x=test_data["resolution"],
        resolution_y=test_data["resolution"],
        nb_thread=1
    )

    # Test X axis
    actual_x = coverage.read_axis_x()

    np.testing.assert_allclose(
        actual_x,
        test_data["target_x"]
    )

    # Test Y axis
    actual_y = coverage.read_axis_y()

    np.testing.assert_allclose(
        actual_y,
        test_data["target_y"]
    )

    # Test bathymetry
    actual_data = coverage.read_variable_bathymetry()

    np.testing.assert_allclose(
        actual_data,
        test_data["target_data"],
        rtol=1e-10,
        atol=1e-10
    )


@pytest.mark.parametrize("interpolator", INTERPOLATORS)
@pytest.mark.parametrize("threads", THREADS)
@pytest.mark.parametrize(
    "test_data",
    [
        (size, resolution, False)  # Test with nonlinear data compare to mono thread
        for size in MATRIX_SIZES
        for resolution in RESOLUTIONS
    ],
    indirect=True
)
def test_interpolation_thread_consistency(
        test_data,
        interpolator,
        threads,
        monkeypatch,
        caplog
):
    caplog.set_level(logging.INFO)

    monkeypatch.setattr(
        Coverage,
        "HORIZONTAL_INTERPOLATOR",
        interpolator
    )

    def calculate(nb_thread):
        reader = MemoryReader(
            x_axis=test_data["source_x"],
            y_axis=test_data["source_y"],
            bathymetry=test_data["data"]
        )

        coverage = Coverage(
            reader=reader,
            resolution_x=test_data["resolution"],
            resolution_y=test_data["resolution"],
            nb_thread=nb_thread
        )

        return coverage.read_variable_bathymetry()

    reference = calculate(1)
    actual = calculate(threads)

    np.testing.assert_allclose(
        actual,
        reference,
        rtol=1e-10,
        atol=1e-10,
        equal_nan=True
    )


@pytest.mark.parametrize("interpolator", ["cgal"])
@pytest.mark.parametrize("threads", THREADS)
@pytest.mark.parametrize("overlap", OVERLAPS)
@pytest.mark.parametrize(
    "test_data",
    [
        (size, resolution, True)  # Test with affine data
        for size in MATRIX_SIZES
        for resolution in RESOLUTIONS
    ],
    indirect=True
)
def test_cgal_multithreaded_interpolation_with_overlap(
        test_data,
        interpolator,
        threads,
        overlap,
        monkeypatch,
        caplog
):
    caplog.set_level(logging.INFO)

    monkeypatch.setattr(
        Coverage,
        "HORIZONTAL_INTERPOLATOR",
        interpolator
    )

    monkeypatch.setattr(
        Coverage,
        "HORIZONTAL_OVERLAPING_SIZE",
        overlap
    )

    def calculate(nb_thread):
        reader = MemoryReader(
            x_axis=test_data["source_x"],
            y_axis=test_data["source_y"],
            bathymetry=test_data["data"]
        )

        coverage = Coverage(
            reader=reader,
            resolution_x=test_data["resolution"],
            resolution_y=test_data["resolution"],
            nb_thread=nb_thread
        )

        return coverage.read_variable_bathymetry()

    actual = calculate(threads)

    np.testing.assert_allclose(
        actual,
        test_data["target_data"],
        rtol=1e-10,
        atol=1e-10,
        equal_nan=True,
        err_msg=(
            f"Interpolation mismatch: "
            f"interpolator={interpolator}, "
            f"threads={threads}, "
            f"overlap={overlap}, "
            f"size={test_data['data'].shape}, "
            f"resolution={test_data['resolution']}"
        )
    )
