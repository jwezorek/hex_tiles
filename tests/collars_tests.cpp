#include "collars.hpp"
#include "collars_detail.hpp"
#include "finite_patch.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
namespace h=jw::hexgrid;namespace g=jw::grid;
void require(bool b,const char* s) {if(!b) throw std::runtime_error(s);}
template<class F> void fails(F f) {try {f();} catch(const std::exception&) {return;}throw std::runtime_error("expected rejection");}
h::collar semantic(const h::hex_patch& p) {h::collar c;for(auto [at,x]:p.tiles)c.push_back({at,x.type,x.orientation});return c;}
h::collar transformed(h::collar c,unsigned turns,h::coordinate translation={},bool mirror=false) {
 for(auto& x:c) {
  const auto m=x.type==g::tile_type::yellow || x.type==g::tile_type::red?2u:6u;
  if(mirror) {x.at={-x.at.q-x.at.r,x.at.r};x.orientation=(m+(m==2?1:3)-x.orientation)%m;}
  for(unsigned t=0;t<turns;++t) x.at={-x.at.r,x.at.q+x.at.r};
  x.orientation=(x.orientation+turns)%m;x.at=x.at+translation;
 }return c;
}
int main() {
 try {
  // Exact byte-for-byte semantic promotion, including rooted context coordinates.
  std::ifstream data(std::string(COLLAR_SOURCE_DIR)+"/research/collar-contexts.txt");
  require(bool(data),"research context fixture");std::string line;size_t vi=0,ci=0;
  while(std::getline(data,line)) {
   if(!line.empty() && line.back()=='\r') line.pop_back();
   if(line.empty() || line[0]=='#') continue;
   std::istringstream row(line);std::string field;std::getline(row,field,';');bool vertex=field=="vertex";
   h::hex_patch p;
   while(std::getline(row,field,';')) {std::replace(field.begin(),field.end(),',',' ');std::istringstream in(field);int q,r,t;unsigned o;in>>q>>r>>t>>o;require(bool(in),"context row");p.insert({{q,r},g::tile_type(t),o});}
   const auto& expected=vertex?h::detail::certified_vertex_seeds().at(vi++):h::detail::certified_radius2_contexts().at(ci++);
   require(semantic(p)==semantic(expected),"compiled certificate differs from research");
  }
  require(vi==19 && ci==29,"certificate cardinalities");
  std::set<h::collar> edges;
  for(const auto& seed:h::detail::certified_vertex_seeds()) {
   for(auto a=seed.tiles.begin();a!=seed.tiles.end();++a) for(auto b=std::next(a);b!=seed.tiles.end();++b) {
    auto raw=h::collar{{a->first,a->second.type,a->second.orientation},{b->first,b->second.type,b->second.orientation}};
    // Edges are unrooted: minimize over the two possible roots as well.
    edges.insert(std::min(h::canonicalize_collar(raw,a->first),h::canonicalize_collar(raw,b->first)));
   }
  }
  require(edges.size()==18,"certified edge count");
  const size_t counts[]{5,14,29,50,77,110,149};const unsigned depths[]{0,1,2,2,3,3,3};
  const unsigned cover[]{0,1,2,3,3,4,4};
  std::array<h::collar_catalog,7> catalogs;
  for(unsigned n=0;n<=6;++n) {
   auto c=h::generate_collars(n);catalogs[n]=c;
   require(c.collars.size()==counts[n],"certified collar count");
   require(c.inflation_depth==depths[n],"minimum catalog depth");
   require(!c.coverage_radius,"catalog certificate must not claim geometric coverage");
   require(std::is_sorted(c.collars.begin(),c.collars.end()),"deterministic ordering");
   require(std::adjacent_find(c.collars.begin(),c.collars.end())==c.collars.end(),"duplicate collars");
   for(auto collar:c.collars) {
    require(collar.size()==1+3*n*(n+1),"ball size");std::set<h::coordinate> positions;
    for(auto x:collar) {
     require(std::max({std::abs(x.at.q),std::abs(x.at.r),std::abs(x.at.q+x.at.r)})<=n,"ball support");positions.insert(x.at);
    }
    require(positions.contains({}) && positions.size()==collar.size(),"root and unique positions");
    for(unsigned t=0;t<6;++t) require(h::canonicalize_collar(transformed(collar,t,{13,-8}),{13,-8})==collar,"proper rotation and translation invariance");
   }
   auto sufficient=h::detail::enumerate_collars(n,cover[n]);
   require(c.collars==sufficient.collars,"catalog/coverage set equality");
  }
  require(h::generate_collars(3).collars==catalogs[3].collars,"repeat enumeration ordering");
  bool chiral=false;for(auto c:catalogs[2].collars) if(h::canonicalize_collar(transformed(c,0,{},true))!=c) chiral=true;
  require(chiral,"reflections must remain distinct");
  const g::integer radii[]{0,1,2,4,11};
  for(unsigned k=0;k<=4;++k) {require(h::detail::coverage_radius(k)==radii[k],"exact R(k)");std::cout<<"R("<<k<<")="<<radii[k]<<std::endl;}
  h::collar_options opt;opt.use_coverage_depth=true;
  auto certified=h::generate_collars(6,opt);require(certified.inflation_depth==4 && certified.coverage_radius==11 && certified.collars==catalogs[6].collars,"coverage API policy");
  auto larger=h::generate_collars(7);require(larger.inflation_depth==4 && larger.coverage_radius==11 && !larger.collars.empty(),"larger-radius exact coverage selection");
  fails([]{h::generate_collars(std::numeric_limits<unsigned>::max());});
  fails([]{h::collar_options o;o.max_depth=1;h::generate_collars(6,o);});
  fails([]{h::collar_options o;o.max_catalog_hexes=4;h::generate_collars(0,o);});
  fails([]{h::collar_options o;o.max_catalog_hexes=100;h::generate_collars(6,o);});
  fails([]{h::canonicalize_collar({});});
  fails([]{h::canonicalize_collar({{{},g::tile_type::red,2}});});
  fails([]{h::canonicalize_collar({{{},g::tile_type::red,0},{{},g::tile_type::red,0}});});
  fails([]{h::canonicalize_collar({{{std::numeric_limits<g::integer>::max(),0},g::tile_type::red,0}},{-1,0});});
  std::cout<<"Catalog counts, support, invariance, certificates and resource errors passed\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
