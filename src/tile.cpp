#include "tile.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {

    jw::line_segment inset_line_segment(const jw::line_segment& seg) {
        // Relative to the arrow segment's length, so even tiny positive
        // custom lengths have a valid inset. Preserve its direction.
        const jw::point delta = seg.v - seg.u;
        const double length = std::hypot(delta.x, delta.y);
        if (!(length > 0.0) || !std::isfinite(length)) {
            throw std::runtime_error("to_svg: arrow segment has zero length");
        }
        const jw::point along = 0.12 * delta;
        const jw::point normal = {-0.06 * delta.y, 0.06 * delta.x};
        const jw::line_segment inset{seg.u + along + normal, seg.v - along + normal};
        if (!std::isfinite(inset.u.x) || !std::isfinite(inset.u.y)
            || !std::isfinite(inset.v.x) || !std::isfinite(inset.v.y)) {
            throw std::overflow_error("to_svg: nonfinite inset arrow");
        }
        return inset;
    }

    void scale_points(std::vector<jw::point>& pts, double scale) {
        for (jw::point& pt : pts) {
            pt.x *= scale;
            pt.y *= scale;
        }
    }

    void scale_points(jw::line_segment& seg, double scale) {
        seg.u *= scale;
        seg.v *= scale;
    }

    std::string svg_escape(const std::string& s) {
        std::string result;
        result.reserve(s.size());

        for (char ch : s) {
            switch (ch) {
            case '&':
                result += "&amp;";
                break;

            case '<':
                result += "&lt;";
                break;

            case '>':
                result += "&gt;";
                break;

            case '"':
                result += "&quot;";
                break;

            case '\'':
                result += "&apos;";
                break;

            default:
                result += ch;
                break;
            }
        }

        return result;
    }

    std::string points_to_svg_points(const std::vector<jw::point>& pts) {
        std::ostringstream out;
        out << std::setprecision(17);

        for (const jw::point& pt : pts) {
            out << pt.x << ',' << pt.y << ' ';
        }

        return out.str();
    }

    std::vector<jw::point> exterior_edge_points(
            const std::vector<jw::point>& vertices,
            const jw::edge_run& run ) {

        if (vertices.empty() || run.count <= 0) {
            return {};
        }

        const int n = static_cast<int>(vertices.size());

        if (run.first < 0 || run.first >= n) {
            throw std::runtime_error("to_svg: exterior edge run has invalid first edge.");
        }

        if (run.count > n) {
            throw std::runtime_error("to_svg: exterior edge run is longer than polygon boundary.");
        }

        std::vector<jw::point> result;
        result.reserve(static_cast<std::size_t>(run.count) + 1);

        for (int k = 0; k <= run.count; ++k) {
            const int vertex_index = (run.first + k) % n;
            result.push_back(vertices[vertex_index]);
        }

        return result;
    }
}

jw::patch jw::inflate(
        const jw::sub_tile_instance& instance,
        const jw::substitution_rules& rules ) {

    const auto rule_iter = rules.find(instance.type);

    if (rule_iter == rules.end()) {
        throw std::runtime_error(
            "inflate: no substitution rule found for subtile type: " + instance.type
        );
    }

    const sub_tile_rule& rule = rule_iter->second;

    patch result;
    result.reserve(rule.inflate_substitution.size());

    for (const sub_tile_instance& child : rule.inflate_substitution) {
        result.push_back({
            .type = child.type,

            // child.trans maps child canonical coordinates into the
            // inflated canonical parent. instance.trans maps the parent
            // into the containing patch. Their composition maps the child
            // directly into the containing patch.
            .trans = instance.trans * child.trans
            });
    }

    return result;
}

jw::patch jw::inflate( const patch& p, const substitution_rules& rules
) {
    patch result;

    std::size_t result_size = 0;

    for (const sub_tile_instance& instance : p) {
        const auto rule_iter = rules.find(instance.type);

        if (rule_iter == rules.end()) {
            throw std::runtime_error(
                "inflate: no substitution rule found for subtile type: " + instance.type
            );
        }

        const std::size_t n = rule_iter->second.inflate_substitution.size();
        if (n > result.max_size() - result_size) {
            throw std::length_error("inflate: patch size overflow");
        }
        result_size += n;
    }

    result.reserve(result_size);

    for (const sub_tile_instance& instance : p) {
        patch inflated_instance = inflate(instance, rules);

        result.insert(
            result.end(),
            std::make_move_iterator(inflated_instance.begin()),
            std::make_move_iterator(inflated_instance.end())
        );
    }

    return result;
}

jw::patch jw::apply_matrix(const matrix& m, const patch& p) {

    patch result;
    result.reserve(p.size());

    for (const sub_tile_instance& instance : p) {
        result.push_back({
            .type = instance.type,
            .trans = m * instance.trans
            });
    }

    return result;
}

std::vector<jw::point> jw::vertices_of(
        const sub_tile_instance& instance, const sub_tile_types& types) {

    const auto type_iter = types.find(instance.type);

    if (type_iter == types.end()) {
        throw std::runtime_error(
            "vertices_of: no subtile type found for instance type: " + instance.type
        );
    }

    const sub_tile_type& type = type_iter->second;

    return apply_matrix(instance.trans, type.vertices);
}

std::vector<jw::tile> jw::to_full_tiles( const patch& p, const sub_tile_types& types ) {
    static_cast<void>(p);
    static_cast<void>(types);

    // TODO:
    // Reconstruct full tiles by grouping subtile instances whose internal
    // edges match and whose external edge runs form the boundary of a full
    // tile.
    return {};
}

void jw::to_svg(const std::string& out_path, const patch& tiles,
        const sub_tile_types& tile_dict,
        bool show_arrows,
        double scale, double margin) {

    if (!(scale > 0.0) || !std::isfinite(scale)) {
        throw std::runtime_error("to_svg: scale must be greater than zero.");
    }

    if (!(margin >= 0.0) || !std::isfinite(margin)) {
        throw std::runtime_error("to_svg: margin must be non-negative.");
    }

    struct drawable_sub_tile {
        std::string type;
        std::string color;
        std::vector<point> vertices;
        std::vector<point> exterior;
        line_segment arrow{};
    };

    std::vector<drawable_sub_tile> drawable_tiles;
    drawable_tiles.reserve(tiles.size());

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (const sub_tile_instance& instance : tiles) {
        const auto type_iter = tile_dict.find(instance.type);

        if (type_iter == tile_dict.end()) {
            throw std::runtime_error(
                "to_svg: no subtile type found for instance type: " + instance.type
            );
        }

        const sub_tile_type& type = type_iter->second;
        std::vector<point> vertices = apply_matrix(instance.trans, type.vertices);
        scale_points(vertices, scale);
        line_segment arrow{};
        if (show_arrows) {
            if (type.arrow.u < 0 || type.arrow.v < 0
                || static_cast<std::size_t>(type.arrow.u) >= type.vertices.size()
                || static_cast<std::size_t>(type.arrow.v) >= type.vertices.size()) {
                throw std::runtime_error("to_svg: invalid arrow indices for " + instance.type);
            }
            arrow = inset_line_segment({type.vertices[type.arrow.u],
                                        type.vertices[type.arrow.v]});
            arrow = apply_matrix(instance.trans, arrow);
            scale_points(arrow, scale);
            if (!std::isfinite(arrow.u.x) || !std::isfinite(arrow.u.y)
                || !std::isfinite(arrow.v.x) || !std::isfinite(arrow.v.y)) {
                throw std::overflow_error("to_svg: nonfinite transformed arrow");
            }
        }

        if (vertices.size() < 3) {
            throw std::runtime_error(
                "to_svg: subtile type has fewer than three vertices: " + instance.type
            );
        }

        for (const point& pt : vertices) {
            if (!std::isfinite(pt.x) || !std::isfinite(pt.y)) {
                throw std::overflow_error("to_svg: nonfinite scaled vertex");
            }
            min_x = std::min(min_x, pt.x);
            min_y = std::min(min_y, pt.y);
            max_x = std::max(max_x, pt.x);
            max_y = std::max(max_y, pt.y);
        }

        auto exterior = exterior_edge_points(vertices, type.external_edges);
        

        drawable_tiles.push_back({
            .type = instance.type,
            .color = type.color,
            .vertices = std::move(vertices),
            .exterior = std::move(exterior),
            .arrow = arrow
        });
    }

    if (drawable_tiles.empty()) {
        min_x = 0.0;
        min_y = 0.0;
        max_x = 100.0;
        max_y = 100.0;
    }

    const double width = max_x - min_x;
    const double height = max_y - min_y;
    const double svg_width = width + 2.0 * margin;
    const double svg_height = height + 2.0 * margin;
    if (!(svg_width > 0.0 && svg_height > 0.0)
        || !std::isfinite(svg_width) || !std::isfinite(svg_height)) {
        throw std::overflow_error("to_svg: nonfinite SVG bounds");
    }

    const auto to_svg_point = [&](const point& pt) -> point {
        return {
            (pt.x - min_x) + margin,
            (max_y - pt.y) + margin
        };
        };

    for (drawable_sub_tile& tile : drawable_tiles) {
        for (point& pt : tile.vertices) {
            pt = to_svg_point(pt);
        }

        for (point& pt : tile.exterior) {
            pt = to_svg_point(pt);
        }

        if (show_arrows) {
            tile.arrow.u = to_svg_point(tile.arrow.u);
            tile.arrow.v = to_svg_point(tile.arrow.v);
        }
    }

    std::ofstream out(out_path);

    if (!out) {
        throw std::runtime_error("to_svg: could not open output file: " + out_path);
    }

    out << std::setprecision(17);
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << svg_width << "\" height=\"" << svg_height << "\" "
        << "viewBox=\"0 0 " << svg_width << ' ' << svg_height << "\">\n";

    if (show_arrows) {
        out << "  <defs>\n";
        out << "    <marker id=\"arrowhead\" markerWidth=\"4\" markerHeight=\"3\" "
            << "refX=\"3.8\" refY=\"1.5\" orient=\"auto\" markerUnits=\"strokeWidth\">\n";
        out << "      <polygon points=\"0 0, 4 1.5, 0 3\" fill=\"#000000\" />\n";
        out << "    </marker>\n";
        out << "  </defs>\n";
    }

    out << "  <g id=\"subtiles\">\n";

    for (const drawable_sub_tile& tile : drawable_tiles) {
        out << "    <g class=\"subtile\" data-type=\""
            << svg_escape(tile.type)
            << "\">\n";

        out << "      <polygon points=\""
            << points_to_svg_points(tile.vertices)
            << "\" fill=\""
            << svg_escape(tile.color)
            << "\" stroke=\"none\" />\n";

        if (tile.exterior.size() >= 2) {
            out << "      <polyline points=\""
                << points_to_svg_points(tile.exterior)
                << "\" fill=\"none\" stroke=\"#000000\" "
                << "stroke-width=\"1\" stroke-linecap=\"round\" "
                << "stroke-linejoin=\"round\" />\n";
        }

        if (show_arrows) {
            out << "      <line "
                << "x1=\"" << tile.arrow.u.x << "\" "
                << "y1=\"" << tile.arrow.u.y << "\" "
                << "x2=\"" << tile.arrow.v.x << "\" "
                << "y2=\"" << tile.arrow.v.y << "\" "
                << "stroke=\"#000000\" "
                << "stroke-width=\"1.5\" "
                << "stroke-linecap=\"round\" "
                << "marker-end=\"url(#arrowhead)\" />\n";
        }

        out << "    </g>\n";
    }

    out << "  </g>\n";
    out << "</svg>\n";
    out.close();
    if (!out) {
        throw std::runtime_error("to_svg: writing output failed: " + out_path);
    }
}
