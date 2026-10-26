#include "collars.hpp"
#include "collars_detail.hpp"
#include "finite_patch.hpp"
#include <algorithm>
#include <set>
#include <stdexcept>

namespace jw::hexgrid {
namespace {
    coordinate rotate(coordinate p) { return {grid::subtract(0,p.r),grid::add(p.q,p.r)}; }
    unsigned period(grid::tile_type t) {
        if(static_cast<unsigned>(t)>=grid::types.size()) throw std::invalid_argument("invalid collar color");
        return t==grid::tile_type::yellow || t==grid::tile_type::red ? 2u:6u;
    }
    std::vector<coordinate> ball(unsigned n) {
        if(n>32) throw std::length_error("collar radius exceeds supported limit 32");
        const grid::integer r=n;
        const auto count=grid::add(1,grid::multiply(3,grid::multiply(r,grid::add(r,1))));
        std::vector<coordinate> result;result.reserve(static_cast<size_t>(count));
        for(auto q=-r;q<=r;++q) for(auto s=std::max(-r,-q-r);s<=std::min(r,-q+r);++s) result.push_back({q,s});
        return result;
    }
    // Exact rule counts before substituting any seed. Saturate at the cell cap.
    void preflight(unsigned depth) {
        if(depth>5) throw std::length_error("collar inflation depth exceeds supported limit 5");
        const auto rules=grid::rules({});
        std::array<grid::integer,5> counts{1,3,3,3,1};
        for(unsigned k=0;k<depth;++k) {
            std::array<grid::integer,5> next{};
            for(size_t i=0;i<5;++i) for(auto c:rules[i])
                next[i]=std::min(grid::integer(grid::max_patch_cells)+1,grid::add(next[i],counts[unsigned(c.type)]));
            counts=next;
        }
        for(const auto& seed:detail::certified_vertex_seeds()) {
            grid::integer sum=0;
            for(const auto& [at,x]:seed.tiles) sum=grid::add(sum,grid::multiply(counts[unsigned(x.type)],period(x.type)==2?6:2));
            if(sum>grid::integer(grid::max_patch_cells)) throw std::length_error("collar seed exceeds triangle cell budget");
        }
    }
}
collar canonicalize_collar(const collar& input,coordinate center) {
    if(input.empty()) throw std::invalid_argument("cannot canonicalize an empty collar");
    collar current;current.reserve(input.size());
    for(auto x:input) {
        if(x.orientation>=period(x.type)) throw std::invalid_argument("noncanonical collar orientation");
        x.at=x.at-center;current.push_back(x);
    }
    collar best;
    for(unsigned t=0;t<6;++t) {
        auto key=current;std::sort(key.begin(),key.end());
        for(size_t i=1;i<key.size();++i) if(key[i-1].at==key[i].at) throw std::invalid_argument("duplicate collar position");
        if(t==0 || key<best) best=std::move(key);
        for(auto& x:current) { x.at=rotate(x.at);x.orientation=(x.orientation+1)%period(x.type); }
    }
    return best;
}
namespace detail {
collar_catalog enumerate_collars(unsigned radius,unsigned depth,const collar_options& options) {
    const auto offsets=ball(radius);preflight(depth);
    const auto budget=std::min(options.max_catalog_hexes,size_t(1'000'000));
    if(depth>options.max_depth) throw std::length_error("collar depth exceeds requested budget");
    if(offsets.size()>budget) throw std::length_error("collar storage budget too small");
    std::set<collar> unique;
    for(const auto& seed:certified_vertex_seeds()) {
        // Retain the entire triangular patch between inflations.
        const auto triangles=substitute_triangles(to_triangles(seed),depth);
        const auto hexes=convert(triangles,std::nullopt,reconstruction::supported_sectors);
        for(const auto& [center,tile]:hexes.tiles) {
            collar c;c.reserve(offsets.size());
            for(auto d:offsets) {
                const auto* x=hexes.find(center+d);if(!x) break;
                c.push_back({d,x->type,x->orientation});
            }
            if(c.size()!=offsets.size()) continue;
            auto key=canonicalize_collar(c);
            if(!unique.contains(key)) {
                if(unique.size()>=budget/offsets.size()) throw std::length_error("collar catalog exceeds storage budget");
                unique.insert(std::move(key));
            }
        }
    }
    collar_catalog result;result.radius=radius;result.inflation_depth=depth;
    result.collars.assign(unique.begin(),unique.end());return result;
}
}
collar_catalog generate_collars(unsigned radius,const collar_options& options) {
    const auto offsets=ball(radius);
    if(offsets.size()>std::min(options.max_catalog_hexes,size_t(1'000'000)))
        throw std::length_error("collar storage budget too small");
    if(radius==0) {
        if(options.max_catalog_hexes<5) throw std::length_error("collar storage budget too small");
        collar_catalog result;for(auto t:grid::types) result.collars.push_back({{{},t,0}});return result;
    }
    if(!options.use_coverage_depth && radius<=6) {
        constexpr unsigned depths[]{0,1,2,2,3,3,3};
        auto result=detail::enumerate_collars(radius,depths[radius],options);
        result.depth_basis=collar_depth_basis::certified_catalog;return result;
    }
    // S_k > 4n proves a finite search bound; do not presume monotonic R(k).
    grid::lengths dims;unsigned bound=0;
    while(dims.S<=grid::multiply(4,radius)) {dims=grid::next_lengths(dims);++bound;}
    for(unsigned k=0;k<=std::min({bound,options.max_depth,5u});++k) {
        const auto r=detail::coverage_radius(k);
        if(r>=radius) {
            auto result=detail::enumerate_collars(radius,k,options);
            result.depth_basis=collar_depth_basis::coverage_certificate;result.coverage_radius=r;return result;
        }
    }
    throw std::length_error("coverage certificate requires a depth beyond the resource budget");
}
}
