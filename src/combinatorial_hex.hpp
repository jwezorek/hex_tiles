#pragma once

#include "triangle_grid.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <string_view>
#include <vector>

namespace jw::grid {
    enum class tile_type { yellow, blue, green, purple, red };
    inline constexpr std::array types{tile_type::yellow,tile_type::blue,tile_type::green,
                                      tile_type::purple,tile_type::red};
    std::string_view name(tile_type type);
    std::string_view color(tile_type type);
    tile_type parse_type(std::string_view label);

    struct lengths {
        integer L=2, S=1;
        integer middle() const { return subtract(L,S); }
        auto operator<=>(const lengths&) const = default;
    };
    lengths next_lengths(lengths value);
    struct child { tile_type type; affine transform; };
    using rule_set = std::array<std::vector<child>,5>;
    rule_set rules(lengths value);
    // Exact rows of elementary cells: q>=0, r>=0, q+r<=base, r<=height.
    std::vector<cell> canonical_region(tile_type type, lengths value = {});

    using node_id = std::size_t;
    inline constexpr node_id no_parent = std::numeric_limits<node_id>::max();
    inline constexpr std::size_t max_patch_cells = 1'000'000;
    struct occurrence {
        tile_type type;
        unsigned level;
        node_id parent;
        std::size_t child_slot;
        affine transform; // canonical local frame -> root frame, including reflection
        std::vector<node_id> children;
        std::vector<cell> cells; // empty for non-leaves; owned cells for level zero
    };
    struct incidence { cell triangle; node_id owner; };
    using edge_index = std::map<edge,std::vector<incidence>>;
    struct patch {
        // Stable deterministic preorder IDs, with root ID 0. References into the
        // vector need not survive construction; IDs do. Child order matches hex.cpp.
        std::vector<occurrence> nodes;
        std::vector<node_id> leaves;
        std::map<cell,node_id> occupancy;
        // One incident cell = perimeter; two different owners = subtile seam.
        edge_index edges() const;
    };
    patch make_supertile(tile_type seed, unsigned depth);
}
