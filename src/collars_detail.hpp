#pragma once
// Internal certificate machinery, exposed here for reproducible verification.
#include "collars.hpp"
namespace jw::hexgrid::detail {
    const std::vector<hex_patch>& certified_vertex_seeds();
    const std::vector<hex_patch>& certified_radius2_contexts();
    grid::integer coverage_radius(unsigned depth);
    collar_catalog enumerate_collars(unsigned radius, unsigned depth, const collar_options& options={});
}
