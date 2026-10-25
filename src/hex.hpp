#pragma once

#include "tile.hpp"

#include <cstddef>
#include <string>
#include <string_view>

namespace jw {

    // Leaf dimensions: long trapezoid base L, small length S, and M = L - S.
    // The historical name "small" does not imply S <= M.
    struct hex_lengths {
        double L;
        double S;

        [[nodiscard]] double middle() const { return L - S; }
    };

    [[nodiscard]] hex_lengths golden_lengths();
    [[nodiscard]] hex_lengths regular_lengths();
    void validate_lengths(hex_lengths lengths);
    [[nodiscard]] hex_lengths next_lengths(hex_lengths lengths);

    [[nodiscard]] sub_tile_types get_hex_subtiles(hex_lengths lengths);

    // RAW placements: children at 'lengths' tile the corresponding canonical
    // parent at next_lengths(lengths). There is no inflation scale factor.
    [[nodiscard]] substitution_rules get_hex_rules(hex_lengths lengths);

    // 'lengths' describes the leaves at every depth. Depth zero is one
    // identity instance. Ordered bottom-up composition is child * leaf.
    [[nodiscard]] patch make_supertile(
        std::string_view seed_type, unsigned depth, hex_lengths lengths);

    void rules_svg(const std::string& out_path, hex_lengths lengths,
                   bool show_arrows = true);
}
