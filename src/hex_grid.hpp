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
        // Six CCW sectors beginning at angle -150deg. These auxiliary arrays
        // are derived provenance only; color, orientation, coordinate and
        // patch lattice suffice for reconstructing the underlying subtiles.
        std::array<grid::cell,6> triangles{};
        std::array<grid::node_id,6> source_leaves{};
        // Bit k means sector k was observed during conversion. Unobserved
        // source_leaves are no_parent; triangles in those slots are unspecified.
        unsigned observed_sectors = 0;
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

    // Recognize six-sector motifs from leaf triangle occupancy. A rooted
    // supertile supplies an enclosing expected support; a finite leaf-only
    // patch instead validates the full support of each prototile occurrence.
    // Invalid complete hexes and missing subtile cells are errors; partial
    // boundary hexes are discarded. Fully absent hexes may form intentional
    // holes in a finite patch. Specify placement to preserve the phase.
    // Reconstruct leaf prototiles using only axial positions, types and orientations.
    // Auxiliary triangle/source arrays in hex_tile are not used as inputs.
    grid::patch to_triangles(const hex_patch& hexes);

    // The optional placement specifies the phase exactly, avoiding phase inference
    // on finite patches. Without it, the original supertile behavior is preserved.
    enum class reconstruction { complete_hexes, supported_sectors };
    // supported_sectors decodes every represented hex from one or more sectors.
    // All base prototiles must be complete and agree on phase, including at
    // boundaries. It never infers hexes with no observed sector.
    hex_patch convert(const grid::patch& triangles,
                      std::optional<lattice> placement = std::nullopt,
                      reconstruction mode = reconstruction::complete_hexes);
}
