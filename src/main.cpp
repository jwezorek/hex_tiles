#include "hex.hpp"
#include "geometry.hpp"
#include "combinatorial_hex.hpp"
#include "combinatorial_svg.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
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
            << "  --arrows                 Draw direction arrows (off by default)\n"
            << "  --no-arrows              Disable direction arrows\n"
            << "  --help                   Show this help\n"
            << "Use either --preset or both --long and --short, never both.\n"
            << "Combinatorial mode requires --preset regular; custom lengths are unsupported.\n"
            << "Depth 0 is one canonical leaf; depths >= 1 use raw level-specific rules.\n";
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

    unsigned parse_depth(std::string_view value) {
        unsigned depth = 0;
        const auto* first = value.data();
        const auto [end, error] = std::from_chars(first, first + value.size(), depth);
        if (error != std::errc{} || end != first + value.size()) {
            throw std::invalid_argument("invalid nonnegative integer for --depth: "
                                        + std::string(value));
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
            if (arg == "--arrows") {
                result.arrows = true;
                continue;
            }
            if (arg == "--no-arrows") {
                result.arrows = false;
                continue;
            }
            if (arg != "--preset" && arg != "--long" && arg != "--short"
                && arg != "--seed" && arg != "--depth" && arg != "--out-dir" && arg != "--mode") {
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
                result.depth = parse_depth(value);
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
            std::cout << "Wrote " << patch.leaves.size() << " leaf subtiles ("
                      << patch.occupancy.size() << " triangle cells) to " << opt.out_dir
                      << " (patch.svg, rules.svg)\n";
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
