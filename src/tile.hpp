#pragma once

#include "geometry.hpp"

#include <string>
#include <unordered_map>
#include <vector>

/*------------------------------------------------------------------------------------------------*/

namespace jw {

    struct edge_run {
        int first = 0;
        int count = 0;
    };

    struct line_seg_indices {
        int u;
        int v;
    };

    struct sub_tile_type {
        std::string label;

        // Canonical vertices for this subtile type.
        //
        // Edge i is:
        //
        //     vertices[i] -> vertices[(i + 1) % vertices.size()]
        //
        std::vector<point> vertices;

        // Consecutive canonical boundary edges that lie on the boundary
        // of the full tile containing this subtile.
        //
        // All other edges are internal seams between subtiles of that
        // full tile.
        edge_run external_edges;
        line_seg_indices arrow;
        std::string color;
    };

    using sub_tile_types = std::unordered_map<std::string, sub_tile_type>;

    struct sub_tile_instance {
        // Key into sub_tile_types and substitution_rules.
        std::string type;

        // Transformation from this subtile type's canonical coordinates
        // into the coordinates of the patch containing this instance.
        matrix trans;
    };

    using patch = std::vector<sub_tile_instance>;

    // A reconstructed full tile, represented as a group of subtile instances.
    using tile = std::vector<sub_tile_instance>;

    struct sub_tile_rule {
        // Key into sub_tile_types.
        std::string base_type;

        // Raw child placements from leaf dimensions (L,S) into the parent
        // at next_lengths(L,S). No scaling/normalization is applied.
        // Each transform maps the child's canonical coordinates into the
        // next-level parent's coordinates.
        patch inflate_substitution;
    };

    using substitution_rules = std::unordered_map<std::string, sub_tile_rule>;

    // Low-level ONE-STEP application only. For multiple levels using raw
    // hex rules, use make_supertile: the lengths change at every level and
    // repeatedly applying one fixed rules dictionary is generally invalid.
    patch inflate(
        const sub_tile_instance& instance,
        const substitution_rules& rules
    );

    patch inflate(
        const patch& p,
        const substitution_rules& rules
    );

    patch apply_matrix(const matrix& m, const patch& p);

    std::vector<point> vertices_of(
        const sub_tile_instance& instance,
        const sub_tile_types& types
    );

    std::vector<tile> to_full_tiles(
        const patch& p,
        const sub_tile_types& types
    );

    void to_svg(const std::string& out_path, const patch& tiles,
        const sub_tile_types& tile_dict,
        bool show_arrows,
        double scale, double margin
    );

}