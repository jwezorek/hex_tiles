#include "hex.hpp"
#include "geometry.hpp"
#include "combinatorial_hex.hpp"
#include "combinatorial_svg.hpp"
#include "hex_grid.hpp"
#include "hex_grid_svg.hpp"
#include "finite_patch.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace {
    struct options {
        std::string mode = "geometric";
        jw::hex_lengths lengths = jw::golden_lengths();
        std::string seed = "green";
        unsigned depth = 4;
        std::filesystem::path out_dir = "output";
        bool arrows = false;
        bool hex_svg = false;
        bool finite_patch = false;
        unsigned patch_radius = 2;
        unsigned substitute_steps = 2;
        bool help = false;
    };

    void usage(std::ostream& out) {
        out << "Usage: hex_tiles [options]\n"
            << "  --mode geometric|combinatorial  Implementation (default: geometric)\n"
            << "  --preset golden|regular   Leaf lengths (default: golden)\n"
            << "  --long L --short S       Custom leaf lengths, L > S > 0\n"
            << "  --seed TYPE              yellow|blue|green|purple|red (default: green)\n"
            << "  --depth N                Nonnegative substitution depth (default: 4)\n"
            << "  --out-dir PATH           Directory for patch.svg and rules.svg (default: output)\n"
            << "  --arrows                 Draw subtile direction arrows (off by default)\n"
            << "  --hex-svg                Also write hexagons.svg with orientation arrows\n"
            << "  --finite-patch           Select/interior-substitute a finite hex patch\n"
            << "  --patch-radius N         Radius of selected hex ball (default: 2)\n"
            << "  --substitute-steps N     Steps for selected finite patch (default: 2)\n"
            << "  --no-arrows              Disable direction arrows\n"
            << "  --help                   Show this help\n"
            << "Use either --preset or both --long and --short, never both.\n"
            << "Combinatorial mode requires --preset regular; custom lengths are unsupported.\n"
            << "Depth 0 is one canonical leaf; depths >= 1 use raw level-specific rules.\n"
            << "Finite patch mode writes finite_original.svg and finite_substituted.svg.\n";
    }

    double parse_double(std::string_view value, std::string_view option) {
        double number = 0.0;
        const auto* first = value.data();
        const auto* last = first + value.size();
        const auto [end, error] = std::from_chars(first, last, number,
                                                 std::chars_format::general);
        if (error != std::errc{} || end != last || !std::isfinite(number)) {
            throw std::invalid_argument("invalid finite number for " + std::string(option)
                                        + ": " + std::string(value));
        }
        return number;
    }

    unsigned parse_depth(std::string_view value,std::string_view option="--depth") {
        unsigned depth = 0;
        const auto* first = value.data();
        const auto [end, error] = std::from_chars(first, first + value.size(), depth);
        if (error != std::errc{} || end != first + value.size()) {
            throw std::invalid_argument("invalid nonnegative integer for " + std::string(option)
                                        + ": " + std::string(value));
        }
        return depth;
    }

    options parse(int argc, char* argv[]) {
        options result;
        bool has_preset = false, has_long = false, has_short = false;
        double L = 0.0, S = 0.0;
        std::string preset;
        for (int i = 1; i < argc; ++i) {
            const std::string_view arg(argv[i]);
            if (arg == "--help") {
                result.help = true;
                return result;
            }
            if (arg == "--finite-patch") {
                result.finite_patch = true;
                continue;
            }
            if (arg == "--hex-svg") {
                result.hex_svg = true;
                continue;
            }
            if (arg == "--arrows") {
                result.arrows = true;
                continue;
            }
            if (arg == "--no-arrows") {
                result.arrows = false;
                continue;
            }
            if (arg != "--preset" && arg != "--long" && arg != "--short"
                && arg != "--seed" && arg != "--depth" && arg != "--out-dir" && arg != "--mode"
                && arg != "--patch-radius" && arg != "--substitute-steps") {
                throw std::invalid_argument("unknown option: " + std::string(arg));
            }
            if (++i == argc) {
                throw std::invalid_argument("missing value for " + std::string(arg));
            }
            const std::string_view value(argv[i]);
            if (arg == "--mode") {
                if (value != "geometric" && value != "combinatorial")
                    throw std::invalid_argument("unknown mode: " + std::string(value));
                result.mode = value;
            } else if (arg == "--preset") {
                if (has_preset) throw std::invalid_argument("--preset specified twice");
                has_preset = true;
                preset = value;
            } else if (arg == "--long") {
                if (has_long) throw std::invalid_argument("--long specified twice");
                has_long = true;
                L = parse_double(value, arg);
            } else if (arg == "--short") {
                if (has_short) throw std::invalid_argument("--short specified twice");
                has_short = true;
                S = parse_double(value, arg);
            } else if (arg == "--seed") {
                result.seed = value;
            } else if (arg == "--depth") {
                result.depth = parse_depth(value,arg);
            } else if (arg == "--patch-radius") {
                result.finite_patch = true;
                result.patch_radius = parse_depth(value,arg);
            } else if (arg == "--substitute-steps") {
                result.finite_patch = true;
                result.substitute_steps = parse_depth(value,arg);
            } else {
                result.out_dir = std::filesystem::path(std::string(value));
            }
        }
        if (has_preset && (has_long || has_short)) {
            throw std::invalid_argument("--preset cannot be combined with --long/--short");
        }
        if (has_long != has_short) {
            throw std::invalid_argument("custom dimensions require both --long and --short");
        }
        if (result.hex_svg && result.mode != "combinatorial")
            throw std::invalid_argument("--hex-svg requires --mode combinatorial");
        if (result.finite_patch && result.mode != "combinatorial")
            throw std::invalid_argument("finite patch substitution requires --mode combinatorial");
        if (result.mode == "combinatorial") {
            if (has_long || has_short)
                throw std::invalid_argument("combinatorial mode does not support custom lengths; use --preset regular");
            if (!has_preset || preset != "regular")
                throw std::invalid_argument("combinatorial mode requires --preset regular");
        }
        if (has_preset) {
            if (preset == "golden") result.lengths = jw::golden_lengths();
            else if (preset == "regular") result.lengths = jw::regular_lengths();
            else throw std::invalid_argument("unknown preset: " + preset);
        } else if (has_long) {
            result.lengths = {L, S};
        }
        jw::validate_lengths(result.lengths);
        if (result.seed != "yellow" && result.seed != "blue" && result.seed != "green"
            && result.seed != "purple" && result.seed != "red") {
            throw std::invalid_argument("unknown seed type: " + result.seed);
        }
        if (result.out_dir.empty()) {
            throw std::invalid_argument("--out-dir cannot be empty");
        }
        return result;
    }
    // Find a full radius-r ball, preferring one near the patch centroid.
    // The input supertile is only a convenient DEMO source; the substitute
    // operation below receives the isolated hexes and no ancestry.
    jw::hexgrid::hex_patch select_hex_ball(const jw::hexgrid::hex_patch& source,unsigned radius) {
        namespace h=jw::hexgrid;
        namespace g=jw::grid;
        const auto r=static_cast<g::integer>(radius);
        const auto count=g::add(1,g::multiply(3,g::multiply(r,g::add(r,1))));
        if (count>static_cast<g::integer>(source.size()))
            throw std::invalid_argument("no complete hex ball of requested radius; increase --depth");
        long double mean_q=0,mean_r=0;
        for (const auto& [at,tile]:source.tiles) {
            (void)tile;
            mean_q+=static_cast<long double>(at.q);
            mean_r+=static_cast<long double>(at.r);
        }
        mean_q/=source.size();mean_r/=source.size();
        std::optional<h::coordinate> best;
        long double best_distance=std::numeric_limits<long double>::infinity();
        for (const auto& [at,tile]:source.tiles) {
            (void)tile;
            std::size_t included=0;
            for (const auto& [other,unused]:source.tiles) {
                (void)unused;
                const auto dq=g::subtract(other.q,at.q),dr=g::subtract(other.r,at.r);
                if (dq>=-r && dq<=r && dr>=-r && dr<=r &&
                    g::add(dq,dr)>=-r && g::add(dq,dr)<=r) ++included;
            }
            const bool complete=included==static_cast<std::size_t>(count);
            const long double q=static_cast<long double>(at.q)-mean_q;
            const long double v=static_cast<long double>(at.r)-mean_r;
            const long double score=q*q+q*v+v*v;
            if (complete && score<best_distance) {best=at;best_distance=score;}
        }
        if (!best)
            throw std::invalid_argument("no complete interior hex ball; increase --depth or reduce --patch-radius");
        h::hex_patch selected;
        selected.placement=source.placement;
        for (const auto& [at,tile]:source.tiles) {
            const auto dq=g::subtract(at.q,best->q),dr=g::subtract(at.r,best->r);
            if (dq>=-r && dq<=r && dr>=-r && dr<=r &&
                g::add(dq,dr)>=-r && g::add(dq,dr)<=r) selected.insert(tile);
        }
        return selected;
    }

}

int main(int argc, char* argv[]) {
    try {
        const options opt = parse(argc, argv);
        if (opt.help) {
            usage(std::cout);
            return 0;
        }
        if (opt.mode == "combinatorial") {
            const auto patch = jw::grid::make_supertile(jw::grid::parse_type(opt.seed), opt.depth);
            std::filesystem::create_directories(opt.out_dir);
            jw::grid::to_svg((opt.out_dir / "patch.svg").string(), patch, opt.arrows);
            jw::grid::rules_svg((opt.out_dir / "rules.svg").string(), opt.arrows);
            if (opt.hex_svg || opt.finite_patch) {
                const auto hexes = jw::hexgrid::convert(patch);
                if (opt.hex_svg)
                    jw::hexgrid::to_svg((opt.out_dir / "hexagons.svg").string(), hexes);
                std::cout << "Converted " << hexes.size() << " complete hexagons ("
                          << hexes.rejected_boundary << " boundary fragments discarded)\n";
                if (opt.finite_patch) {
                    const auto selected=select_hex_ball(hexes,opt.patch_radius);
                    const auto result=jw::hexgrid::substitute(selected,opt.substitute_steps);
                    jw::hexgrid::to_svg((opt.out_dir / "finite_original.svg").string(),selected);
                    jw::hexgrid::to_svg((opt.out_dir / "finite_substituted.svg").string(),result);
                    std::cout << "Substituted " << selected.size() << " isolated hexes "
                              << opt.substitute_steps << " steps to " << result.size()
                              << " complete hexes (" << result.rejected_boundary
                              << " boundary fragments discarded)\n";
                }
            }
            std::cout << "Wrote " << patch.leaves.size() << " leaf subtiles ("
                      << patch.occupancy.size() << " triangle cells) to " << opt.out_dir
                      << " (patch.svg, rules.svg"
                      << (opt.hex_svg ? ", hexagons.svg" : "")
                      << (opt.finite_patch ? ", finite_original.svg, finite_substituted.svg" : "")
                      << ")\n";
            return 0;
        }
        const auto tiles = jw::get_hex_subtiles(opt.lengths);
        const auto patch = jw::make_supertile(opt.seed, opt.depth, opt.lengths);

        // The seed's outer canonical polygon at the 'depth'-updated dimensions
        // determines only the SVG display scale, never the mathematical transforms.
        auto outer = opt.lengths;
        for (unsigned level = 0; level < opt.depth; ++level) {
            outer = jw::next_lengths(outer);
        }
        const auto outer_tiles = jw::get_hex_subtiles(outer);
        const auto& parent_vertices = outer_tiles.at(opt.seed).vertices;
        double span = 0.0;
        for (const auto& vertex : parent_vertices) {
            span = std::max(span, std::max(vertex.x, vertex.y));
        }
        if (!(span > 0.0) || !std::isfinite(span)) {
            throw std::overflow_error("nonfinite or zero patch bounds");
        }
        const double display_scale = 1200.0 / span;
        if (!(display_scale > 0.0) || !std::isfinite(display_scale)) {
            throw std::overflow_error("nonfinite SVG display scale");
        }

        std::filesystem::create_directories(opt.out_dir);
        jw::to_svg((opt.out_dir / "patch.svg").string(), patch, tiles,
                   opt.arrows, display_scale, 20.0);
        jw::rules_svg((opt.out_dir / "rules.svg").string(), opt.lengths,
                      opt.arrows);
        std::cout << "Wrote " << patch.size() << " leaf subtiles to "
                  << opt.out_dir << " (patch.svg, rules.svg)\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "hex_tiles: " << error.what() << "\n";
        usage(std::cerr);
        return 1;
    }
}
