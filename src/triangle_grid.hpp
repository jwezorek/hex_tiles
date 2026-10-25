#pragma once

#include <array>
#include <compare>
#include <cstdint>

namespace jw::grid {
    using integer = std::int64_t;
    integer add(integer a, integer b);
    integer subtract(integer a, integer b);
    integer multiply(integer a, integer b);

    // q*(1,0) + r*(1/2,sqrt(3)/2). No Euclidean coordinates in this module.
    struct vertex {
        integer q = 0, r = 0;
        auto operator<=>(const vertex&) const = default;
    };
    vertex operator+(vertex a, vertex b);
    vertex operator-(vertex a, vertex b);

    struct edge {
        vertex a, b; // sorted, so shared edges have the same identity
        edge(vertex u, vertex v);
        auto operator<=>(const edge&) const = default;
    };

    struct cell {
        vertex anchor;
        bool down = false;
        auto operator<=>(const cell&) const = default;
        // Counterclockwise in the lattice basis.
        std::array<vertex, 3> vertices() const;
        std::array<edge, 3> edges() const;
        std::array<cell, 3> neighbors() const;
    };

    // Integer affine maps. a*b means apply b first, then a.
    struct affine {
        integer a = 1, b = 0, c = 0, d = 1;
        vertex translation{};
        auto operator<=>(const affine&) const = default;
        vertex operator()(vertex p) const;
        cell operator()(cell t) const; // rejects maps that do not preserve cells
        integer determinant() const;
        affine inverse() const; // requires determinant +1 or -1
    };
    affine operator*(const affine& a, const affine& b);
    affine translate(vertex offset);
    affine rotate(int sixth_turns);
    affine reflect(); // Euclidean x -> -x, y -> y
}
