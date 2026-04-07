#!/bin/bash
# Download JPL ephemeris files and Earth orientation parameters.
#
# de440s.bsp          (~31 MB)  — required by all programs
# de441_part-1.bsp    (~1.5 GB) — needed for create_de441s (1849–1969)
# de441_part-2.bsp    (~1.5 GB) — needed for create_de441s (1969–2150)
# finals2000A.data    (~2.5 MB) — IERS EOP (polar motion, UT1-UTC)
# latest_eop2.long    (~3.8 MB) — JPL EOP2 (same source as Horizons)

set -euo pipefail
cd "$(dirname "$0")"

NAIF_URL="https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets"
IERS_URL="https://maia.usno.navy.mil/ser7"
JPL_EOP_URL="https://eop2-external.jpl.nasa.gov/eop2"

download() {
    local url="$1"
    local file="$2"
    if [ -f "$file" ]; then
        echo "$file already exists, skipping"
    else
        echo "Downloading $file ..."
        curl -fSL -o "$file.part" "$url/$file"
        mv "$file.part" "$file"
        echo "$file done"
    fi
}

download "$NAIF_URL" de440s.bsp
download "$NAIF_URL" de441_part-1.bsp
download "$NAIF_URL" de441_part-2.bsp
download "$IERS_URL" finals2000A.data
download "$JPL_EOP_URL" latest_eop2.long
