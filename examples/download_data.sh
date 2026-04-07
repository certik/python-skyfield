#!/bin/bash
# Download JPL ephemeris files needed by the Fortran programs.
#
# de440s.bsp       (~31 MB)  — required by all programs
# de441_part-1.bsp (~1.5 GB) — needed for create_de441s (1849–1969)
# de441_part-2.bsp (~1.5 GB) — needed for create_de441s (1969–2150)
#                               and test_de441_horizons

set -euo pipefail
cd "$(dirname "$0")"

NAIF_URL="https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets"

download() {
    local file="$1"
    if [ -f "$file" ]; then
        echo "$file already exists, skipping"
    else
        echo "Downloading $file ..."
        curl -fSL -o "$file.part" "$NAIF_URL/$file"
        mv "$file.part" "$file"
        echo "$file done"
    fi
}

download de440s.bsp
download de441_part-1.bsp
download de441_part-2.bsp
