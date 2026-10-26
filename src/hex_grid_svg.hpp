#pragma once
#include "hex_grid.hpp"
#include "collars.hpp"
#include <string>

namespace jw::hexgrid {
    // Every tile has an outline and an orientation arrow, independent of the
    // old --arrows switch for subtiles. An empty patch produces an empty SVG.
    void to_svg(const std::string& path,const hex_patch& patch);
    struct collar_svg_options {
        std::size_t columns=5;
        double hex_size=12; // circumradius in SVG user units
    };
    // Standalone contact sheet; reads only the supplied catalog. Empty catalogs
    // produce a labeled empty sheet. Invalid input and write failures throw.
    void collars_to_svg(const std::string& path,const collar_catalog& catalog,
                        const collar_svg_options& options={});
}
