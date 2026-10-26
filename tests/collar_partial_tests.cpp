#include "hex_grid.hpp"
#include <algorithm>
#include <bit>
#include <iostream>
#include <stdexcept>
namespace h=jw::hexgrid;namespace g=jw::grid;
void require(bool x,const char* m) {if(!x) throw std::runtime_error(m);}
template<class F> void fails(F f) {try {f();}catch(const std::logic_error&){return;}throw std::runtime_error("expected malformed input rejection");}
h::hex_patch partial(const g::patch& p,std::optional<h::lattice> phase={}) {return h::convert(p,phase,h::reconstruction::supported_sectors);}
int main() {
 try {
  for(auto type:g::types) for(unsigned o=0;o<((type==g::tile_type::yellow || type==g::tile_type::red)?2u:6u);++o) {
   h::hex_patch p;p.placement.origin={8,-7};p.insert({{2,-3},type,o});auto full=h::to_triangles(p);
   for(auto id:full.leaves) {
    g::patch part;part.nodes.push_back(full.nodes[id]);part.leaves.push_back(0);
    for(auto c:part.nodes[0].cells) part.occupancy.emplace(c,0);
    auto out=partial(part,p.placement);
    require(out.size()==1,"a complete base prototile must recover its supported hex");
    const auto& tile=out.at({2,-3});
    require(tile.type==type && tile.orientation==o,"one sector uniquely determines label");
    require(std::popcount(tile.observed_sectors)==part.occupancy.size(),"actual sector provenance");
    for(unsigned i=0;i<6;++i) {
     if(tile.observed_sectors&(1u<<i)) require(part.occupancy.contains(tile.triangles[i]) && tile.source_leaves[i]==0,"observed provenance");
     else require(tile.source_leaves[i]==g::no_parent,"absent provenance must not be fabricated");
    }
    auto inferred=partial(part);require(inferred.size()==1 && inferred.tiles.begin()->second.orientation==o,"frame phase inference");
    require(h::convert(part,p.placement).empty(),"default complete-only behavior retained");
    fails([&]{partial(part,h::lattice{p.placement.origin+g::vertex{1,0}});});
    auto bad=part;bad.nodes[0].transform.a=2;fails([&]{partial(bad);});
    bad=part;bad.occupancy.begin()->second=1000;fails([&]{partial(bad);});
    bad=part;bad.nodes[0].cells.pop_back();fails([&]{partial(bad);});
    bad=part;bad.occupancy.erase(bad.occupancy.begin());fails([&]{partial(bad);});
   }
   auto out=partial(full,p.placement);require(out.size()==1 && out.at({2,-3}).observed_sectors==63,"full observation");
  }
  h::hex_patch p;p.insert({{},g::tile_type::yellow,0});auto tri=h::to_triangles(p);
  // Replace an entire valid one-cell owner with the opposite decoration/color.
  auto incompatible=h::hex_patch{};incompatible.insert({{},g::tile_type::yellow,1});auto other=h::to_triangles(incompatible);
  auto bad=tri;const auto cell=bad.nodes[0].cells[0];bad.nodes[0]=other.nodes[other.occupancy.at(cell)];
  fails([&]{partial(bad);});
  bad=tri;bad.nodes[0].type=g::tile_type::red;fails([&]{partial(bad);});
  // Two individually complete prototiles on incompatible phases.
  bad=tri;auto extra=tri.nodes[0];extra.transform=g::translate({30,0})*extra.transform;
  for(auto& c:extra.cells)c=g::translate({30,0})(c);
  // 31 changes residue; choose it explicitly after the large separating move.
  extra.transform=g::translate({1,0})*extra.transform;for(auto& c:extra.cells)c=g::translate({1,0})(c);
  auto id=bad.nodes.size();bad.nodes.push_back(extra);bad.leaves.push_back(id);for(auto c:extra.cells)bad.occupancy.emplace(c,id);
  fails([&]{partial(bad);});
  h::hex_patch hole;for(auto at:h::directions)hole.insert({at,g::tile_type::yellow,0});
  auto out=partial(h::to_triangles(hole));require(out.size()==6 && !out.find({}),"intentional absent hex must remain absent");
  std::cout<<"Supported sectors, provenance, phase, boundaries and malformed observations passed\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
