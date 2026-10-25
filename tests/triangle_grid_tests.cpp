#include "triangle_grid.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace jw::grid;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void overflow(F action) {
    try { action(); } catch (const std::overflow_error&) { return; }
    throw std::runtime_error("expected integer overflow");
}
int main() {
    try {
        require(rotate(1)(vertex{1,0}) == vertex{0,1}, "60 degree basis rotation");
        require(rotate(1)(vertex{0,1}) == vertex{-1,1}, "60 degree second basis");
        require(reflect()(vertex{2,3}) == vertex{-5,3}, "vertical reflection");
        require(rotate(6) == affine{}, "six rotations");
        require(reflect()*reflect() == affine{}, "two reflections");
        for (int k = -6; k <= 6; ++k) {
            for (bool mirror : {false, true}) {
                const auto t = translate({-7,11}) * rotate(k) * (mirror ? reflect() : affine{});
                require(t.inverse()*t == affine{}, "left inverse");
                require(t*t.inverse() == affine{}, "right inverse");
                require(reflect()*rotate(k)*reflect() == rotate(-k), "dihedral relation");
                for (integer q = -2; q <= 2; ++q) for (integer r = -2; r <= 2; ++r) {
                    for (bool down : {false, true}) {
                        const cell x{{q,r},down};
                        require(t.inverse()(t(x)) == x, "cell round trip");
                        for (const auto n : x.neighbors()) {
                            auto a = t(x).edges(), b = t(n).edges();
                            std::ptrdiff_t shared = 0;
                            for (auto e : a) shared += std::count(b.begin(),b.end(),e);
                            require(shared == 1, "adjacency preserved");
                        }
                    }
                }
                for (int j=0; j<6; ++j) {
                    const vertex v{-3,5};
                    require((t*rotate(j))(v) == t(rotate(j)(v)), "composition order");
                }
            }
        }
        constexpr auto hi = std::numeric_limits<integer>::max();
        constexpr auto lo = std::numeric_limits<integer>::min();
        overflow([&]{ (void)add(hi,1); });
        overflow([&]{ (void)subtract(lo,1); });
        overflow([&]{ (void)multiply(lo,-1); });
        overflow([&]{ (void)multiply(hi,2); });
        overflow([&]{ (void)translate({hi,0})(vertex{1,0}); });
        require(multiply(lo,1)==lo && multiply(hi,-1)==-hi, "valid extreme products");
        require(subtract(lo,lo)==0 && add(lo,hi)==-1, "valid extreme sums");
        bool rejected = false;
        try { (void)affine{2,0,0,2,{}}(cell{}); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "non-cell-preserving map rejected");
        std::cout << "Triangle lattice transformations passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
