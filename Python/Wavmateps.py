# -*- coding: utf-8 -*-
"""Wavmateps.py
Python translation of the MATLAB function `Wavmat_eps.m`.

Constructs the N×N orthogonal transform matrix for an ε-decimated
orthogonal discrete wavelet transform (periodic boundary conditions).

- h      : orthonormal low-pass filter (row or column)
- N      : signal length (power of 2)
- eps    : length-k0 vector with entries in {0,1}; per-level phase choice
- shift  : integer alignment parameter (default 2)

Notes
-----
- If eps = ones(k0), this matches the "standard" phase convention used in
  the original Wavmat layout.
- The inverse transform is W.T (transpose), since W is orthonormal.
"""

"""
Wavmateps.py


Example
-------
>>> import numpy as np
>>> dat = np.array([1, 0, -3, 2, 1, 0, 1, 2], dtype=float)
>>> h   = np.array([np.sqrt(2)/2, np.sqrt(2)/2])
>>> W1  = Wavmat_eps(h, 8, [1, 1, 1], shift=2)
>>> wt1 = W1 @ dat
>>> rec1 = W1.T @ wt1
>>> np.max(np.abs(rec1 - dat))
0.0
"""

from __future__ import annotations

from typing import Iterable, Optional, Union

import numpy as np


ArrayLike = Union[np.ndarray, Iterable[float]]


def Wavmateps(h: ArrayLike, N: int, eps: ArrayLike, shift: int = 2) -> np.ndarray:
    """
    Build the N×N orthogonal transform matrix for an ε-decimated orthogonal DWT (periodic).

    Parameters
    ----------
    h : array-like
        Orthonormal low-pass filter coefficients (1D).
    N : int
        Signal length (must be a power of 2).
    eps : array-like
        Length-k0 vector with entries in {0,1}. Each entry selects the decimation phase
        at the corresponding level:
            eps[k]=1 -> "odd" phase
            eps[k]=0 -> "even" phase (phase-flipped decimation)
    shift : int, default 2
        Integer alignment parameter .

    Returns
    -------
    W : (N, N) ndarray
        Orthogonal transform matrix; inverse is W.T.

    Raises
    ------
    ValueError
        If N is not a power of 2, or if eps has invalid length/values.
    """
    if not isinstance(N, (int, np.integer)) or N <= 0:
        raise ValueError("N must be a positive integer.")

    # Check N is power of 2
    if N & (N - 1) != 0:
        raise ValueError("N must be a power of 2 (e.g., 256, 1024, ...).")

    # Levels
    J = int(np.log2(N))

    # eps handling
    eps_arr = np.asarray(list(eps), dtype=int).ravel()
    if eps_arr.size == 0:
        raise ValueError("eps must be a non-empty 1D array-like of 0/1 values.")
    eps_arr = np.mod(eps_arr, 2)
    k0 = int(eps_arr.size)
    if k0 > J:
        raise ValueError(f"eps length k0={k0} cannot exceed J=log2(N)={J}.")

    # low-pass filter
    h0 = np.asarray(h, dtype=complex).ravel()

    # Make QMF high-pass from low-pass
    signs = (-1.0) ** np.arange(1, h0.size + 1)
    g0 = np.flip(np.conj(h0) * signs)

    # Zero-extend to length N for modular indexing
    if h0.size < N:
        hpad = np.concatenate([h0, np.zeros(N - h0.size, dtype=complex)])
    else:
        hpad = h0[:N].copy()

    if g0.size < N:
        gpad = np.concatenate([g0, np.zeros(N - g0.size, dtype=complex)])
    else:
        gpad = g0[:N].copy()

    # Initialization: identity in the coarsest space
    oldmat = np.eye(2 ** (J - k0), dtype=complex)

    # Descend levels: k = k0, k0-1, ..., 1
    for k in range(k0, 0, -1):
        ubJk = 2 ** (J - k)       # number of columns in the new h/g blocks
        ubJk1 = 2 ** (J - k + 1)  # number of rows in the new h/g blocks

        # Per-level phase control via an effective shift:
        # shift_k = shift + (1 - eps(k));
        shift_k = int(shift + (1 - eps_arr[k - 1]))

        # Build hmat, gmat (with per-level shift_k)
        hmat = np.zeros((ubJk1, ubJk), dtype=complex)
        gmat = np.zeros((ubJk1, ubJk), dtype=complex)


        for jj0 in range(ubJk):        # jj0 = 0..ubJk-1  corresponds to jj = jj0+1
            for ii0 in range(ubJk1):   # ii0 = 0..ubJk1-1 corresponds to ii = ii0+1
                idx = (N + (ii0 + 1) - 2 * (jj0 + 1) + shift_k - 1) % ubJk1
                hmat[ii0, jj0] = hpad[idx]
                gmat[ii0, jj0] = gpad[idx]

        # Assemble next-level transform:
        top = oldmat @ hmat.T
        bottom = gmat.T
        W = np.vstack([top, bottom])

        oldmat = W

    # If the input is real, return a real matrix (up to numerical noise)
    if np.all(np.isreal(h0)) and np.all(np.isreal(np.asarray(eps_arr))):
        W = np.real_if_close(W, tol=1000)

    return W

import numpy as np
dat = np.array([1, 0, -3, 2, 1, 0, 1, 2], dtype=float)
h   = np.array([np.sqrt(2)/2, np.sqrt(2)/2])
W1  = Wavmat_eps(h, 8, [1, 1, 1], shift=2)
wt1 = W1 @ dat
rec1 = W1.T @ wt1
np.max(np.abs(rec1 - dat))
#0.0
