#include "hex.hpp"
#include "geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
    constexpr double h = 0.86602540378443864676372317075293618;
    constexpr double sixty_degrees = std::numbers::pi / 3.0;
    constexpr jw::point O{0.0, 0.0};
    constexpr std::array<const char*, 5> labels = {
        "yellow", "blue", "green", "purple", "red"
    };
    constexpr std::size_t max_cached_leaves = 1'000'000;

    double checked(double value, const char* what) {
        if (!std::isfinite(value)) {
            throw std::overflow_error(std::string("nonfinite geometry: ") + what);
        }
        return value;
    }

    jw::point pt(double x, double y) {
        return {checked(x, "x coordinate"), checked(y, "y coordinate")};
    }

    jw::sub_tile_instance place(const char* type, jw::point from, jw::point to,
                                int k, bool mirror = false) {
        using namespace jw;
        // 'from' is specified in the reflected coordinate system.
        const matrix transform = map_pt_to_pt_and_rotate(from, to, k * sixty_degrees)
            * (mirror ? scale_matrix(-1.0, 1.0) : identity_matrix());
        if (!transform.allFinite()) {
            throw std::overflow_error("nonfinite child placement transform");
        }
        return {type, transform};
    }

    struct bounds {
        double min_x = std::numeric_limits<double>::infinity();
        double min_y = std::numeric_limits<double>::infinity();
        double max_x = -std::numeric_limits<double>::infinity();
        double max_y = -std::numeric_limits<double>::infinity();

        void add(jw::point p) {
            checked(p.x, "transformed vertex x");
            checked(p.y, "transformed vertex y");
            min_x = std::min(min_x, p.x);
            min_y = std::min(min_y, p.y);
            max_x = std::max(max_x, p.x);
            max_y = std::max(max_y, p.y);
        }
        double width() const { return checked(max_x - min_x, "bounds width"); }
        double height() const { return checked(max_y - min_y, "bounds height"); }
    };

    bounds patch_bounds(const jw::patch& patch, const jw::sub_tile_types& tiles) {
        bounds b;
        for (const auto& instance : patch) {
            for (auto vertex : jw::vertices_of(instance, tiles)) {
                b.add(vertex);
            }
        }
        return b;
    }
}

jw::hex_lengths jw::golden_lengths() {
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    return {phi, 1.0 / phi};
}

jw::hex_lengths jw::regular_lengths() {
    return {2.0, 1.0};
}

void jw::validate_lengths(hex_lengths lengths) {
    if (!std::isfinite(lengths.L) || !std::isfinite(lengths.S)
        || !(lengths.L > lengths.S && lengths.S > 0.0)) {
        throw std::invalid_argument("lengths must be finite and satisfy L > S > 0");
    }
    const double M = lengths.middle();
    if (!(M > 0.0) || !std::isfinite(M)) {
        throw std::invalid_argument("the derived middle length must be positive and finite");
    }
    checked(h * lengths.L, "leaf height");
    checked(h * lengths.S, "leaf height");
    checked(h * M, "leaf height");
}

jw::hex_lengths jw::next_lengths(hex_lengths lengths) {
    validate_lengths(lengths);
    // Use sums of positive quantities rather than 3*L-S, which can
    // overflow needlessly. A = 2L-S = L+M; D = 3L-S = L+A.
    const double A = checked(lengths.L + lengths.middle(), "next middle length");
    const double D = checked(lengths.L + A, "next long base");
    checked(h * A, "next parent height");
    if (!(D > lengths.L && A > 0.0)) {
        throw std::overflow_error("substitution lengths lost precision");
    }
    return {D, lengths.L};
}

jw::sub_tile_types jw::get_hex_subtiles(hex_lengths lengths) {
    validate_lengths(lengths);
    const double L = lengths.L, S = lengths.S, M = lengths.middle();
    std::vector<sub_tile_type> tiles = {
        {"yellow", equilateral_triangle_vertices(M), {1, 1}, {0, 2}, "#ffff00"},
        {"blue", equilateral_trapezoid_vertices(L, M), {0, 3}, {0, 3}, "#3fbcec"},
        {"green", equilateral_trapezoid_vertices(L, M), {0, 3}, {0, 3}, "#3da347"},
        {"purple", equilateral_trapezoid_vertices(L, S), {0, 3}, {0, 3}, "#6327ac"},
        {"red", equilateral_triangle_vertices(S), {1, 1}, {0, 2}, "#ff0000"}
    };
    sub_tile_types result;
    for (auto& tile : tiles) {
        const auto label = tile.label;
        result.emplace(label, std::move(tile));
    }
    return result;
}

jw::substitution_rules jw::get_hex_rules(hex_lengths lengths) {
    // Also checks all parent dimensions before constructing placements.
    const hex_lengths parent = next_lengths(lengths);
    const double L = lengths.L, M = lengths.middle();
    const double A = parent.middle(), D = parent.L;
    const point yellow_C = pt(M + L / 2.0, checked(h * L, "yellow C"));
    const point blue_C = pt(M / 2.0, checked(h * M, "blue C"));
    const point blue_Q = pt(L + M / 2.0, checked(h * M, "blue Q"));
    const point blue_E = pt(L + M, 0.0);
    const point blue_U = pt(D - A / 2.0, checked(h * A, "blue U"));
    const point blue_V = pt(D, 0.0);
    const point green_C1 = pt(M, 0.0);
    const point green_U = pt(A / 2.0, checked(h * A, "green U"));
    const point green_P = pt(A, 0.0);
    const point green_C2 = pt(D - L / 2.0, checked(h * L, "green C2"));
    const point purple_P = pt(L, 0.0);
    const point purple_U = pt(D - L / 2.0, checked(h * L, "purple U"));
    const point purple_C = pt(purple_U.x - M, purple_U.y);
    const point purple_V = pt(D, 0.0);
    const point red_T = pt(L / 2.0, checked(h * L, "red T"));

    substitution_rules result;
    result.emplace("yellow", sub_tile_rule{"yellow", {
        place("blue",   {L, 0.0}, O, -2),
        place("purple", O, {M, 0.0}, 0),
        place("yellow", O, yellow_C, 2),
        place("yellow", O, yellow_C, -2),
        place("yellow", O, yellow_C, 1, true)
    }});
    result.emplace("blue", sub_tile_rule{"blue", {
        place("yellow", O, blue_C, 0),
        place("yellow", O, blue_C, -2),
        place("yellow", O, blue_C, 3, true),
        place("blue", O, {M, 0.0}, 0),
        place("purple", O, blue_Q, 2),
        place("purple", O, blue_Q, 5, true),
        place("blue", {-L, 0.0}, blue_E, 1, true),
        place("blue", {L, 0.0}, blue_E, -2),
        place("yellow", O, blue_U, -2),
        place("yellow", O, blue_U, 1, true),
        place("red", O, blue_V, 2)
    }});
    result.emplace("green", sub_tile_rule{"green", {
        place("yellow", O, green_C1, 2),
        place("yellow", O, green_C1, -1, true),
        place("yellow", O, green_C1, 0),
        place("green", O, green_U, -2),
        place("green", O, green_U, 3, true),
        place("purple", O, green_P, 0),
        place("purple", O, green_P, 2),
        place("purple", O, green_P, -1, true),
        place("yellow", O, green_C2, 2),
        place("yellow", O, green_C2, -2),
        place("yellow", O, green_C2, 1, true)
    }});
    result.emplace("purple", sub_tile_rule{"purple", {
        place("yellow", O, O, 0),
        place("purple", O, purple_P, 0),
        place("purple", O, purple_P, 2),
        place("purple", O, purple_P, -1, true),
        place("yellow", O, purple_C, -2),
        place("yellow", O, purple_C, 1, true),
        place("yellow", O, purple_C, 3, true),
        place("green", O, purple_V, 2)
    }});
    result.emplace("red", sub_tile_rule{"red", {
        place("yellow", O, red_T, -2),
        place("purple", O, O, 0)
    }});
    return result;
}

jw::patch jw::make_supertile(std::string_view seed_type, unsigned depth,
                             hex_lengths lengths) {
    validate_lengths(lengths);
    const std::string seed(seed_type);
    if (!std::any_of(labels.begin(), labels.end(),
                     [&](const char* name) { return seed_type == name; })) {
        throw std::invalid_argument("unknown seed type: " + seed);
    }

    std::unordered_map<std::string, patch> previous;
    for (const char* label : labels) {
        previous.emplace(label, patch{{label, identity_matrix()}});
    }

    for (unsigned level = 0; level < depth; ++level) {
        // The patch at this level has canonical geometry at 'lengths';
        // never reuse a previous level's raw rules.
        const auto rules = get_hex_rules(lengths);
        std::unordered_map<std::string, patch> next;
        std::size_t total = 0;
        for (const char* parent : labels) {
            const auto& children = rules.at(parent).inflate_substitution;
            std::size_t count = 0;
            for (const auto& child : children) {
                const std::size_t n = previous.at(child.type).size();
                if (n > max_cached_leaves - total - count) {
                    throw std::length_error("requested depth exceeds the 1,000,000 cached-leaf limit");
                }
                count += n;
            }
            total += count;
            patch composed;
            composed.reserve(count);
            for (const auto& child : children) {
                for (const auto& leaf : previous.at(child.type)) {
                    const matrix transform = child.trans * leaf.trans;
                    if (!transform.allFinite()) {
                        throw std::overflow_error("nonfinite composed leaf transform");
                    }
                    composed.push_back({leaf.type, transform});
                }
            }
            next.emplace(parent, std::move(composed));
        }
        previous = std::move(next);
        lengths = next_lengths(lengths);
    }
    return std::move(previous.at(seed));
}

void jw::rules_svg(const std::string& out_path, hex_lengths lengths,
                   bool show_arrows) {
    const auto tiles = get_hex_subtiles(lengths);
    const auto rules = get_hex_rules(lengths);
    patch diagram;
    const double column_gap = checked(0.6 * lengths.L, "rules column gap");
    const double row_gap = checked(0.45 * lengths.L, "rules row gap");
    double cursor = 0.0;

    for (const char* label : labels) {
        const patch one_leaf{{label, identity_matrix()}};
        const patch& children = rules.at(label).inflate_substitution;
        const bounds leaf_b = patch_bounds(one_leaf, tiles);
        const bounds child_b = patch_bounds(children, tiles);
        const double row_height = std::max(leaf_b.height(), child_b.height());
        const double child_x = checked(leaf_b.width() + column_gap, "rules column x");

        // Coordinate y increases upward; negative cursor stacks rows top-down
        // after SVG's usual y-axis inversion.
        const auto append_column = [&](const patch& source, const bounds& b,
                                       double left_x) {
            const double dx = checked(left_x - b.min_x, "rules column placement x");
            const double dy = checked(-cursor - (row_height - b.height()) / 2.0
                                      - b.max_y, "rules column placement y");
            for (const auto& instance : source) {
                const matrix transform = translation_matrix(dx, dy) * instance.trans;
                if (!transform.allFinite()) {
                    throw std::overflow_error("nonfinite rules layout transform");
                }
                diagram.push_back({instance.type, transform});
            }
        };
        append_column(one_leaf, leaf_b, 0.0);
        append_column(children, child_b, child_x);
        cursor = checked(cursor + row_height + row_gap, "rules row placement");
    }
    const bounds diagram_b = patch_bounds(diagram, tiles);
    const double scale = 1300.0 / std::max(diagram_b.width(), diagram_b.height());
    to_svg(out_path, diagram, tiles, show_arrows, checked(scale, "rules display scale"), 20.0);
}
