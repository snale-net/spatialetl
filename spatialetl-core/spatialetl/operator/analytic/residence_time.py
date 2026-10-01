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


def compute_residence_time(sourceAxis, data):
    """
    """

    # --- Dimensions ---
    Nt, Ny, Nx = data.shape
    Npix = Nx * Ny

    # --- Reshape to (Nt, Npix) ---
    C = data.reshape(Nt, Npix)

    # --- Initialize ---
    TR_vec = np.full(Npix, np.nan)

    # --- Mask valid pixels ---
    valid = ~np.all(np.isnan(C), axis=0)
    C_valid = C[:, valid]

    if C_valid.size == 0:
        return TR_vec.reshape(Ny, Nx)

    # --- Threshold ---
    seuil = C_valid[0, :] / np.exp(1)

    # --- Find crossings ---
    below = C_valid <= seuil[None, :]  # (Nt, Nvalid)

    # First True index along time axis
    idx = np.argmax(below, axis=0)

    # Detect pixels that never cross
    never_cross = ~np.any(below, axis=0)

    # Valid crossings
    valid_cross = (idx > 0) & (~never_cross)

    # --- Prepare output for valid pixels ---
    TR_valid = np.full(C_valid.shape[1], np.nan)

    # --- Interpolation ---
    i = np.where(valid_cross)[0]
    idx_i = idx[i]

    c0 = C_valid[idx_i - 1, i]
    c1 = C_valid[idx_i, i]
    t0 = sourceAxis[idx_i - 1]
    t1 = sourceAxis[idx_i]

    # Linear interpolation
    TR_valid[i] = t0 + (seuil[i] - c0) * (t1 - t0) / (c1 - c0)

    # --- Edge cases ---
    TR_valid[~valid_cross] = sourceAxis[-1]

    # --- Put back ---
    TR_vec[valid] = TR_valid

    # --- Reshape ---
    return TR_vec.reshape(Ny, Nx)


def compute_lagrangian_residence_time(sourceAxis, data):
    TR_lag = np.full(data.shape[2], np.max(sourceAxis))

    dt = np.median(np.diff(sourceAxis))
    t0 = sourceAxis[0]

    # Search indices where tracer vanish (=> lon, lat is Nan)
    mask = np.isnan(data)
    has_nan = mask.any(axis=1)
    idx = np.argmax(has_nan, axis=0)  # First index where Nan are found

    # If no NaN were found, set -1
    no_nan = ~has_nan.any(axis=0)
    idx[no_nan] = -1

    valid = idx != -1  # Create a mask where values are found

    # Select indices where values are found
    tracer_idx = np.nonzero(valid)[0]
    time_idx = idx[valid]

    # Compute the TR
    TR_lag[tracer_idx] = sourceAxis[time_idx] - t0 + dt

    return TR_lag
