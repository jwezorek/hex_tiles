#pragma once
#include "hex_grid.hpp"
#include <string>

namespace jw::hexgrid {
    // Every tile has an outline and an orientation arrow, independent of the
    // old --arrows switch for subtiles. An empty patch produces an empty SVG.
    void to_svg(const std::string& path,const hex_patch& patch);
}
