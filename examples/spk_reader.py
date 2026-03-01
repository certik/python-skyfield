#!/usr/bin/env python3
"""
Pure Python/NumPy implementation of JPL SPK (Spacecraft and Planet Kernel)
file reader.

Implements the NAIF DAF (Double precision Array File) binary format and
SPK Type 2/3 Chebyshev polynomial interpolation.  Produces results that
are bit-for-bit identical to jplephem.

References:
  - https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html
  - https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html
"""

import struct
import numpy as np
from numpy import array, rollaxis

T0 = 2451545.0        # J2000.0 epoch (Julian Date)
S_PER_DAY = 86400.0   # seconds per day


# ═══════════════════════════════════════════════════════════════════════
#  DAF — Double precision Array File reader
# ═══════════════════════════════════════════════════════════════════════

class DAF:
    """Reader for NASA SPICE Double precision Array Files (DAF).

    A DAF file consists of 1024-byte records.  Record 1 is the file
    record (header), followed by optional comment records, then one or
    more pairs of summary/name records.

    Each summary record contains a control area (3 doubles: next, prev,
    count) followed by packed segment descriptors.  The paired name
    record stores the human-readable name of each segment.

    SPK files use nd=2 (two doubles: start/end seconds from J2000)
    and ni=6 (six ints: target, center, frame, data_type, start_i, end_i).
    """

    def __init__(self, file_object):
        self.file = file_object

        file_record = self._read_record(1)

        # Detect endianness by trying both byte orders
        locidw = file_record[:8].upper().rstrip()
        if locidw != b'NAIF/DAF' and not locidw.startswith(b'DAF/'):
            raise ValueError(
                f'Not a DAF file: starts with {locidw!r}')

        for locfmt, endian in ((b'LTL-IEEE', '<'), (b'BIG-IEEE', '>')):
            s = struct.Struct(endian + '8sII60sIII8s603s28s297s')
            values = s.unpack(file_record)
            if values[1] == 2:   # nd == 2 for SPK files
                self.endian = endian
                self.nd = values[1]
                self.ni = values[2]
                self.fward = values[4]   # first summary record
                self.bward = values[5]   # last summary record
                self.free = values[6]    # first free word
                break
        else:
            raise ValueError(
                'neither endianness scan produces the expected ND=2')

        # Build struct formats for summary parsing
        summary_format = 'd' * self.nd + 'i' * self.ni
        self._control_struct = struct.Struct(self.endian + 'ddd')
        self._summary_struct = struct.Struct(self.endian + summary_format)
        length = self._summary_struct.size
        self._summary_step = length + (-length % 8)  # pad to 8 bytes

    # ── low-level I/O ──

    def _read_record(self, n):
        """Read 1024-byte record *n* (1-indexed)."""
        self.file.seek((n - 1) * 1024)
        return self.file.read(1024)

    def read_array(self, start, end):
        """Return doubles from word *start* to *end* inclusive (1-indexed)."""
        length = 1 + end - start
        self.file.seek(8 * (start - 1))
        data = self.file.read(8 * length)
        return np.ndarray(length, self.endian + 'd', data)

    def map_array(self, start, end):
        """Return doubles from word *start* to *end* inclusive (1-indexed).

        Identical to read_array — jplephem uses memory-mapping here, but
        a plain read is simpler and works everywhere.
        """
        return self.read_array(start, end)

    # ── segment enumeration ──

    def summaries(self):
        """Yield ``(name_bytes, descriptor_tuple)`` for every segment."""
        record_number = self.fward
        while record_number:
            data = self._read_record(record_number)
            next_number, _prev, n_summaries = \
                self._control_struct.unpack(data[:24])

            name_data = self._read_record(record_number + 1)

            step = self._summary_step
            ctrl_size = self._control_struct.size

            for i in range(int(n_summaries)):
                j = ctrl_size + i * step
                summary_bytes = data[j:j + self._summary_struct.size]
                descriptor = self._summary_struct.unpack(summary_bytes)

                name = name_data[i * step : (i + 1) * step].strip()
                yield name, descriptor

            record_number = int(next_number)


# ═══════════════════════════════════════════════════════════════════════
#  SPK — Spacecraft and Planet Kernel
# ═══════════════════════════════════════════════════════════════════════

class SPK:
    """A JPL SPK ephemeris kernel for computing positions and velocities.

    Usage::

        kernel = SPK.open('de440s.bsp')
        pos = kernel[0, 3].compute(2451545.0, 0.5)
        pos, vel = kernel[0, 3].compute_and_differentiate(2451545.0, 0.5)
        kernel.close()
    """

    def __init__(self, daf):
        self.daf = daf
        self.segments = [
            Segment(daf, source, descriptor)
            for source, descriptor in daf.summaries()
        ]
        self.pairs = {(s.center, s.target): s for s in self.segments}

    @classmethod
    def open(cls, path):
        """Open an SPK file and return an SPK instance."""
        f = open(path, 'rb')
        try:
            return cls(DAF(f))
        except Exception:
            f.close()
            raise

    def close(self):
        """Close this SPK file."""
        self.daf.file.close()

    def __getitem__(self, key):
        return self.pairs[key]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __str__(self):
        lines = [f'SPK with {len(self.segments)} segments:']
        for s in self.segments:
            lines.append(f'  {s.center} -> {s.target}  Type {s.data_type}')
        return '\n'.join(lines)


# ═══════════════════════════════════════════════════════════════════════
#  Segment — Type 2/3 Chebyshev polynomial interpolation
# ═══════════════════════════════════════════════════════════════════════

class Segment:
    """A single SPK segment backed by Chebyshev polynomials (Type 2 or 3).

    Matches the interface and numerical behaviour of ``jplephem.spk.Segment``
    exactly — including the Clenshaw recurrence, two-part JD handling, and
    the derivative chain-rule factors.
    """

    _data = None   # lazy-loaded

    def __init__(self, daf, source, descriptor):
        self.daf = daf
        self.source = source
        (self.start_second, self.end_second,
         self.target, self.center, self.frame,
         self.data_type, self.start_i, self.end_i) = descriptor
        self.start_jd = T0 + self.start_second / S_PER_DAY
        self.end_jd = T0 + self.end_second / S_PER_DAY

    # ── coefficient loading (mirrors jplephem's Segment._data) ──

    def _load_data(self):
        if self._data is not None:
            return self._data

        if self.data_type == 2:
            component_count = 3
        elif self.data_type == 3:
            component_count = 6
        else:
            raise ValueError(
                f'only SPK types 2 and 3 are supported, got {self.data_type}')

        # Last 4 words of segment: init, intlen, rsize, n
        meta = self.daf.read_array(self.end_i - 3, self.end_i)
        init = meta[0]              # initial epoch, seconds from J2000
        intlen = meta[1]            # interval length, seconds
        rsize = int(meta[2])        # doubles per interval record
        n = int(meta[3])            # number of intervals

        coefficient_count = (rsize - 2) // component_count

        # Read coefficient block (everything except the 4 metadata words)
        coefficients = self.daf.map_array(self.start_i, self.end_i - 4)

        # (n, rsize) → drop MID & RADIUS columns → reshape per component
        coefficients = coefficients.reshape((n, rsize))
        coefficients = coefficients[:, 2:]
        coefficients = coefficients.reshape((n, component_count,
                                             coefficient_count))
        # → (component_count, n, coefficient_count)
        coefficients = rollaxis(coefficients, 1)
        # → (coefficient_count, component_count, n)
        coefficients = rollaxis(coefficients, 2)
        # reverse so highest-degree coefficient is first (Clenshaw order)
        coefficients = coefficients[::-1]

        self._data = (init, intlen, coefficients)
        return self._data

    # ── Clenshaw recurrence (position + differentiation) ──

    def generate(self, tdb, tdb2):
        """Yield *components* then *rates* for time ``tdb + tdb2``.

        This is a generator identical in semantics to
        ``jplephem.spk.Segment.generate``.
        """
        scalar = (not getattr(tdb, 'shape', 0)
                  and not getattr(tdb2, 'shape', 0))
        if scalar:
            tdb = array((tdb,))

        init, intlen, coefficients = self._load_data()
        coefficient_count, component_count, n = coefficients.shape

        # Two-part JD → interval index + offset (seconds within interval)
        index1, offset1 = divmod((tdb - T0) * S_PER_DAY - init, intlen)
        index2, offset2 = divmod(tdb2 * S_PER_DAY, intlen)
        index3, offset = divmod(offset1 + offset2, intlen)
        index = (index1 + index2 + index3).astype(int)

        if (index < 0).any() or (index > n).any():
            raise ValueError('requested time is outside segment coverage')

        # endpoint wrap
        omegas = (index == n)
        index[omegas] -= 1
        offset[omegas] += intlen

        coefficients = coefficients[:, :, index]

        # ── Chebyshev evaluation (Clenshaw) ──
        s = 2.0 * offset / intlen - 1.0
        s2 = 2.0 * s

        w0 = w1 = 0.0
        wlist = []

        for coefficient in coefficients[:-1]:
            w2 = w1
            w1 = w0
            w0 = coefficient + (s2 * w1 - w2)
            wlist.append(w1)

        components = coefficients[-1] + (s * w0 - w1)

        if scalar:
            components = components[:, 0]

        yield components

        # ── Chebyshev differentiation ──
        dw0 = dw1 = 0.0

        for coefficient, w1 in zip(coefficients[:-1], wlist):
            dw2 = dw1
            dw1 = dw0
            dw0 = 2.0 * w1 + dw1 * s2 - dw2

        rates = w0 + s * dw0 - dw1
        rates /= intlen
        rates *= 2.0
        rates *= S_PER_DAY

        if scalar:
            rates = rates[:, 0]

        yield rates

    # ── public API ──

    def compute(self, tdb, tdb2=0.0):
        """Return ``[x, y, z]`` position (km) at ``tdb + tdb2``."""
        for position in self.generate(tdb, tdb2):
            return position

    def compute_and_differentiate(self, tdb, tdb2=0.0):
        """Return ``(position, velocity)`` at ``tdb + tdb2``."""
        return tuple(self.generate(tdb, tdb2))
