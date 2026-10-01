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


def moving_average(x, w):
    """Simple moving average (1D) with window size 2*w+1."""
    kernel = np.ones(2 * w + 1)
    return np.convolve(x, kernel, mode='same')


def nanmoving_average(x, w, interp_nan=False):
    """NaN-aware 1D moving average."""
    x = np.asarray(x, dtype=float)

    valid = np.isfinite(x).astype(float)
    x_filled = np.nan_to_num(x, nan=0.0)

    kernel = np.ones(2 * w + 1)

    sum_x = np.convolve(x_filled, kernel, mode='same')
    count = np.convolve(valid, kernel, mode='same')

    with np.errstate(invalid='ignore', divide='ignore'):
        y = sum_x / count

    if not interp_nan:
        y[~np.isfinite(x)] = np.nan

    return y


def nanmoving_average2(X, m=1, n=None, interp_nan=False):
    """
    2D NaN-aware moving average (MATLAB nanmoving_average2 equivalent).
    Originally written in MATLAB by
    M.S. Carlos Adrián Vargas Aguilera
    Physical Oceanography PhD candidate
    CICESE
    Mexico, october 2006
    nubeobscura@hotmail.com

    Parameters
    ----------
    X : 2D array
    m : window radius (rows)
    n : window radius (cols), if None -> n = m
    interp_nan : bool, if True smooth NaNs too

    Returns
    -------
    Y : smoothed array
    """

    if X.ndim != 2:
        raise ValueError("Input must be a 2D matrix")

    if n is None:
        n = m

    M, N = X.shape
    Y = np.zeros((M, N))

    # Matrix with 1 where valid, 0 where NaN
    A = (~np.isnan(X)).astype(float)

    # --- 1. Column-wise pass ---
    B = np.zeros((M, N))

    for j in range(N):
        Y[:, j] = nanmoving_average(X[:, j], n, interp_nan=True)
        B[:, j] = moving_average(A[:, j], n)

    Y = Y * B  # convert to sums

    # --- 2. Row-wise pass ---
    C = np.zeros((M, N))

    for i in range(M):
        Y[i, :] = nanmoving_average(Y[i, :], m, interp_nan=True)
        B[i, :] = moving_average(B[i, :], m)
        C[i, :] = moving_average(A[i, :], m)

    Y = Y * C  # still sums

    B[B == 0] = 1
    Y = Y / B  # final averaging

    if not interp_nan:
        Y[np.isnan(X)] = np.nan

    return Y
