#pragma once
#include "hex_grid.hpp"
#include <vector>

namespace jw::hexgrid {
    // Semantic data only: no triangle ownership or reconstruction provenance.
    struct collar_hex {
        coordinate at;
        grid::tile_type type;
        unsigned orientation=0;
        auto operator<=>(const collar_hex&) const = default;
    };
    using collar = std::vector<collar_hex>;
    enum class collar_depth_basis { single_colors, certified_catalog, coverage_certificate };
    struct collar_catalog {
        unsigned radius=0;
        unsigned inflation_depth=0;
        collar_depth_basis depth_basis=collar_depth_basis::single_colors;
        // Present only when the selected depth was checked for geometric coverage.
        std::optional<grid::integer> coverage_radius;
        std::vector<collar> collars;
    };
    struct collar_options {
        // Bounded exact implementation: radius <=32, inflation depth <=5,
        // <=1,000,000 triangle cells per seed and <=1,000,000 stored hex labels.
        // Lower budgets may be requested. Exceeding any budget throws length_error.
        unsigned max_depth=5;
        std::size_t max_catalog_hexes=1'000'000;
        // Also compute an independent geometric certificate for n<=6.
        bool use_coverage_depth=false;
    };
    // Translate center to zero and minimize over the six proper rotations.
    // Input orientations must already be canonical (yellow/red <2; others <6).
    // Reflections are distinct. Empty input and duplicate positions are errors.
    collar canonicalize_collar(const collar& input, coordinate center={});
    // Exhaustive legal catalog, ordered by canonical semantic keys. n=0 has
    // five colors. n=1..6 uses certified minimum catalog depths unless requested
    // otherwise; larger radii use the exact 29-context coverage certificate.
    collar_catalog generate_collars(unsigned radius, const collar_options& options={});
}
