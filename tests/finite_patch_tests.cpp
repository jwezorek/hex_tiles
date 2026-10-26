#include "finite_patch.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace g=jw::grid;
namespace h=jw::hexgrid;
namespace {
    void check(bool ok,const char* message) {
        if (!ok) throw std::runtime_error(message);
    }
    template<typename F> void invalid(F f) {
        try {f();}
        catch(const std::invalid_argument&) {return;}
        throw std::runtime_error("invalid input should have been rejected");
    }
    template<typename F> void malformed(F f) {
        try {f();}
        catch(const std::logic_error&) {return;}
        throw std::runtime_error("malformed geometry should have been rejected");
    }
    void equal_hexes(const h::hex_patch& a,const h::hex_patch& b) {
        check(a.placement.origin==b.placement.origin && a.size()==b.size(),"hex placement or count changed");
        auto x=a.tiles.begin(),y=b.tiles.begin();
        for (;x!=a.tiles.end();++x,++y)
            check(x->first==y->first && x->second.type==y->second.type &&
                  x->second.orientation==y->second.orientation,"hex reconstruction changed identity");
    }
    void check_triangle_patch(const g::patch& p) {
        std::set<g::cell> all;
        std::size_t sum=0;
        for (auto id:p.leaves) {
            const auto& n=p.nodes.at(id);
            check(n.level==0 && n.transform.determinant()!=0,"valid leaf affine");
            std::set<g::cell> actual(n.cells.begin(),n.cells.end());
            check(actual.size()==n.cells.size(),"no duplicated cells inside leaf");
            std::set<g::cell> expected;
            for (auto c:g::canonical_region(n.type)) expected.insert(n.transform(c));
            check(actual==expected,"leaf support matches transformed canonical region");
            for (auto c:n.cells) {
                check(p.occupancy.at(c)==id,"triangle occupancy owner");
                check(all.insert(c).second,"triangle overlap");
                ++sum;
            }
        }
        check(sum==p.occupancy.size(),"no missing occupancy or duplicated triangles");
    }
    using key=std::tuple<g::tile_type,g::affine,std::vector<g::cell>>;
    key leaf_key(const g::occurrence& o) {
        auto cells=o.cells;
        std::sort(cells.begin(),cells.end());
        return {o.type,o.transform,std::move(cells)};
    }
    void verify_original_leaves(const g::patch& original,const h::hex_patch& hexes) {
        auto reconstructed=h::to_triangles(hexes);
        check_triangle_patch(reconstructed);
        std::set<g::node_id> leaves;
        for (const auto& [at,tile]:hexes.tiles) {
            (void)at;
            for (auto id:tile.source_leaves) leaves.insert(id);
        }
        std::multiset<key> want,got;
        for (auto id:leaves) want.insert(leaf_key(original.nodes.at(id)));
        for (auto id:reconstructed.leaves) got.insert(leaf_key(reconstructed.nodes.at(id)));
        check(want==got,"reconstructed prototiles differ from original source instances");
    }
    bool matches_region(const h::hex_patch& region,const h::hex_patch& full) {
        if (region.empty()) return true;
        const auto& first=region.tiles.begin()->second;
        for (const auto& [at,candidate]:full.tiles) {
            if (candidate.type!=first.type || candidate.orientation!=first.orientation) continue;
            const auto offset=at-region.tiles.begin()->first;
            bool match=true;
            for (const auto& [loc,tile]:region.tiles) {
                const auto* found=full.find(loc+offset);
                if (!found || found->type!=tile.type || found->orientation!=tile.orientation) {
                    match=false;
                    break;
                }
            }
            if (match) return true;
        }
        return false;
    }
    h::hex_patch ball(const h::hex_patch& src,int radius) {
        h::hex_patch out;
        out.placement=src.placement;
        // For a radius-1 ball, pick the first seven-hex interior neighborhood.
        for (const auto& [at,tile]:src.tiles) {
            (void)tile;
            std::set<h::coordinate> coords{at};
            for (auto n:h::neighbors(at)) coords.insert(n);
            if (radius!=1) throw std::runtime_error("test supports radius 1 only");
            if (std::all_of(coords.begin(),coords.end(),[&](auto c){return src.find(c)!=nullptr;})) {
                for (auto c:coords) out.insert(src.at(c));
                return out;
            }
        }
        throw std::runtime_error("no radius-one ball found");
    }
}
int main() {
    try {
        // Every canonical orientation, including the two parities of each
        // 3-fold-symmetric rosette, must reconstruct without provenance.
        for (auto type:g::types) {
            const unsigned rotations=(type==g::tile_type::yellow || type==g::tile_type::red)?2:6;
            for (unsigned orientation=0;orientation<rotations;++orientation) {
                h::hex_patch one;
                one.placement.origin={-4,7};
                check(one.insert({{2,-3},type,orientation}),"isolated hex insertion");
                auto p=h::to_triangles(one);
                check_triangle_patch(p);
                check(p.leaves.size()==(rotations==2?6u:2u),"isolated hex subtile count");
                equal_hexes(one,h::convert(p,one.placement));
                equal_hexes(one,h::substitute(one,0));
                for (unsigned steps:{1u,2u}) {
                    auto next=h::substitute_triangles(p,steps);
                    check_triangle_patch(next);
                    auto result=h::substitute(one,steps);
                    check(!result.empty(),"nonempty substituted single hex");
                    equal_hexes(result,h::substitute(one,steps));
                    check_triangle_patch(h::to_triangles(result));
                }
            }
        }
        // Global, level-dependent placement is tested against completely
        // independently generated supertiles, with exact coordinate
        // translation registration and disregarding clipped outer rings.
        for (auto type:g::types) {
            const auto original=g::make_supertile(type,2);
            const auto source=h::convert(original);
            verify_original_leaves(original,source);
            equal_hexes(source,h::convert(h::to_triangles(source),source.placement));
            for (unsigned steps:{1u,2u}) {
                const auto expanded=h::substitute(source,steps);
                const auto full=h::convert(g::make_supertile(type,2+steps));
                check(matches_region(expanded,full),"finite substitution disagrees with independent supertile region");
                check_triangle_patch(h::to_triangles(expanded));
            }
        }
        const auto parent=g::make_supertile(g::tile_type::green,3);
        const auto all=h::convert(parent);
        const auto neighborhood=ball(all,1);
        check(neighborhood.size()==7,"connected radius-one neighborhood");
        auto p=h::to_triangles(neighborhood);
        check_triangle_patch(p);
        equal_hexes(neighborhood,h::convert(p,neighborhood.placement));
        for (unsigned steps:{1u,2u}) {
            auto r=h::substitute(neighborhood,steps);
            const auto full=h::convert(g::make_supertile(g::tile_type::green,3+steps));
            check(matches_region(r,full),"mixed hex neighborhood maps into full supertile");
        }
        // Map ordering is intrinsic; insertion order and unused provenance
        // arrays must not influence reconstruction or substitution.
        h::hex_patch reordered;
        reordered.placement=neighborhood.placement;
        for (auto it=neighborhood.tiles.rbegin();it!=neighborhood.tiles.rend();++it) {
            auto t=it->second;
            t.source_leaves.fill(g::no_parent);
            t.triangles.fill({});
            reordered.insert(t);
        }
        equal_hexes(h::substitute(reordered,2),h::substitute(neighborhood,2));
        h::hex_patch empty;
        empty.placement.origin={3,-2};
        equal_hexes(empty,h::substitute(empty,0));
        equal_hexes(empty,h::substitute(empty,2));
        check(h::to_triangles(empty).occupancy.empty(),"empty triangle patch");
        invalid([&] {auto bad=neighborhood;bad.tiles.begin()->second.orientation=7;h::to_triangles(bad);});
        invalid([&] {auto bad=neighborhood;bad.tiles.begin()->second.at={999,888};h::to_triangles(bad);});
        invalid([&] {auto bad=neighborhood;bad.tiles.begin()->second.type=static_cast<g::tile_type>(90);h::to_triangles(bad);});
        invalid([&] {auto bad=p;bad.nodes.front().transform.translation={500,0};h::substitute_triangles(bad,1);});
        invalid([&] {auto bad=p;bad.occupancy.erase(bad.occupancy.begin());h::substitute_triangles(bad,1);});
        // Only an exterior fragment may be dropped. Missing elementary
        // triangles in the middle of an otherwise full seven-hex patch are errors.
        auto damaged=p;
        h::coordinate center{};
        for (const auto& [at,tile]:neighborhood.tiles) {
            (void)tile;
            const auto ns=h::neighbors(at);
            if (std::all_of(ns.begin(),ns.end(),
                [&](auto n){return neighborhood.find(n)!=nullptr;})) {center=at;break;}
        }
        // Avoid iterators from different temporary arrays.
        const auto& central=neighborhood.at(center);
        damaged.occupancy.erase(central.triangles[0]);
        malformed([&]{h::convert(damaged,neighborhood.placement);});
        auto isolated=h::hex_patch{};
        isolated.placement=neighborhood.placement;
        isolated.insert(central);
        // Full substituted leaf occurrences can meet the outer boundary
        // mid-hex; those partial rings must be discarded, without clipping
        // individual source prototiles.
        const auto substituted_triangles=h::substitute_triangles(h::to_triangles(isolated),1);
        const auto substituted_hexes=h::substitute(isolated,1);
        const auto boundary=h::convert(substituted_triangles,substituted_hexes.placement);
        check(boundary.rejected_boundary>0,"incomplete boundary hex discarded");
        equal_hexes(boundary,substituted_hexes);
        // A complete tile removed intentionally from the middle creates
        // an internal patch boundary, not a malformed subtile occurrence.
        h::hex_patch ring;
        ring.placement=neighborhood.placement;
        for (auto n:h::neighbors(center)) ring.insert(neighborhood.at(n));
        equal_hexes(ring,h::convert(h::to_triangles(ring),ring.placement));
        check(!h::substitute(ring,2).empty(),"hex patch with an intentional hole");
        // An intentional disconnected input is rejected because its relative
        // positions after changing the length parameters are not constrained.
        h::hex_patch separated;
        separated.placement=neighborhood.placement;
        auto first=central,second=central;
        second.at=second.at+h::coordinate{40,40};
        separated.insert(first);separated.insert(second);
        invalid([&]{(void)h::substitute(separated,1);});
        std::cout<<"Finite hex reconstruction, orientations, exact substitution, "
                   "boundary checks, invalid inputs and deterministic order passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
