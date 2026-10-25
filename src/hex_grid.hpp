#pragma once

#include "combinatorial_hex.hpp"

#include <array>
#include <compare>
#include <map>
#include <optional>
#include <utility>

namespace jw::hexgrid {
    // Axial basis in the triangle lattice: U=(1,1), V=(-1,2).
    // Hexagon centers are triangle-grid vertices with (q-r) mod 3 fixed.
    // The coordinate type is deliberately distinct from grid::vertex.
    struct coordinate {
        grid::integer q = 0, r = 0;
        auto operator<=>(const coordinate&) const = default;
    };
    coordinate operator+(coordinate a, coordinate b);
    coordinate operator-(coordinate a, coordinate b);
    inline constexpr std::array<coordinate,6> directions{{
        {1,0}, {0,1}, {-1,1}, {-1,0}, {0,-1}, {1,-1}
    }};
    std::array<coordinate,6> neighbors(coordinate at);

    struct lattice {
        // The exact center of axial (0,0), in triangle-grid coordinates.
        grid::vertex origin{};
        grid::vertex center(coordinate at) const;
        coordinate locate(grid::vertex center) const; // rejects off-lattice points
    };

    struct hex_tile {
        coordinate at;
        grid::tile_type type;
        // 0 points along +q (right in Euclidean space); increments rotate 60deg CCW.
        // Yellow/red rosettes have 3-fold symmetry: only orientation parity is distinct.
        // Reflections of all five known motifs are rotation-equivalent; no mirror bit.
        unsigned orientation = 0;
        // Six CCW sectors beginning at angle -150deg, source triangle and leaf ID.
        std::array<grid::cell,6> triangles{};
        std::array<grid::node_id,6> source_leaves{};
    };

    struct hex_patch {
        lattice placement;
        std::map<coordinate,hex_tile> tiles;
        std::size_t rejected_boundary = 0; // incomplete six-cell groups at exterior

        bool insert(hex_tile tile); // false on duplicate coordinate
        const hex_tile* find(coordinate at) const;
        const hex_tile& at(coordinate where) const;
        std::size_t size() const { return tiles.size(); }
        bool empty() const { return tiles.empty(); }
        // Inclusive axial bounds; nullopt for an empty patch.
        std::optional<std::pair<coordinate,coordinate>> bounds() const;
    };

    // Derive a hex patch from leaf triangle occupancy (never used during inflation).
    // Throws on mixed/invalid complete hexes, duplicate triangles, or interior holes.
    // Incomplete hexes on the outside of the root supertile are discarded.
    hex_patch convert(const grid::patch& triangles);
}
