#pragma once

#include "hex_grid.hpp"

namespace jw::hexgrid {
    // Every occurrence in a finite triangle patch is a BASE-LENGTH prototile.
    // Requires compatible shared-corner constraints across its connected
    // components. Returns another base-length leaf patch; no supertile needed.
    grid::patch substitute_triangles(const grid::patch& input, unsigned steps);

    // Full finite-patch pipeline. Only complete boundary hexes are retained.
    hex_patch substitute(const hex_patch& input, unsigned steps);
}
