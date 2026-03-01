#!/usr/bin/env python3
"""Convert nutation.npz to a flat binary file readable by Fortran.

Binary layout (all little-endian, sequential):
  Header: 10 x int32  — array dimensions (rows, cols for 2D; n, 0 for 1D)
  Then each array written as contiguous doubles or int32s.

Order of arrays:
  1. nals_t            678 x 5   int32
  2. lunisolar_lon     678 x 3   float64
  3. lunisolar_obl     678 x 3   float64
  4. napl_t            687 x 14  int32
  5. nut_lon           687 x 2   float64
  6. nut_obl           687 x 2   float64
  7. ke0_t              33 x 14  int32
  8. ke1                14       int32
  9. se0_t_0            33       float64
 10. se0_t_1            33       float64
"""

import numpy as np
import struct
import os

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        '..', 'skyfield', 'data')
nut = np.load(os.path.join(data_dir, 'nutation.npz'))

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nutation.dat')

with open(out, 'wb') as f:
    # Helper: write 2D int array as int32 (Fortran column-major order)
    def write_int2d(arr):
        r, c = arr.shape
        f.write(struct.pack('<ii', r, c))
        f.write(np.asfortranarray(arr.astype('<i4')).tobytes(order='F'))

    # Helper: write 2D float array as float64 (Fortran column-major order)
    def write_float2d(arr):
        r, c = arr.shape
        f.write(struct.pack('<ii', r, c))
        f.write(np.asfortranarray(arr.astype('<f8')).tobytes(order='F'))

    # Helper: write 1D int array as int32
    def write_int1d(arr):
        n = arr.shape[0]
        f.write(struct.pack('<ii', n, 0))
        f.write(arr.astype('<i4').tobytes())

    # Helper: write 1D float array as float64
    def write_float1d(arr):
        n = arr.shape[0]
        f.write(struct.pack('<ii', n, 0))
        f.write(arr.astype('<f8').tobytes())

    write_int2d(nut['nals_t'])                           # 678 x 5
    write_float2d(nut['lunisolar_longitude_coefficients']) # 678 x 3
    write_float2d(nut['lunisolar_obliquity_coefficients']) # 678 x 3
    write_int2d(nut['napl_t'])                            # 687 x 14
    write_float2d(nut['nutation_coefficients_longitude'])  # 687 x 2
    write_float2d(nut['nutation_coefficients_obliquity'])  # 687 x 2
    write_int2d(nut['ke0_t'])                             # 33 x 14
    write_int1d(nut['ke1'])                               # 14
    write_float1d(nut['se0_t_0'])                         # 33
    write_float1d(nut['se0_t_1'])                         # 33

print(f'Wrote {os.path.getsize(out)} bytes to {out}')
print('Arrays:')
for name in ['nals_t', 'lunisolar_longitude_coefficients',
             'lunisolar_obliquity_coefficients', 'napl_t',
             'nutation_coefficients_longitude', 'nutation_coefficients_obliquity',
             'ke0_t', 'ke1', 'se0_t_0', 'se0_t_1']:
    a = nut[name]
    print(f'  {name:45s} shape={str(a.shape):15s} dtype={a.dtype}')
