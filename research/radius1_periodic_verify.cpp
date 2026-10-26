// Independent C++ check of a periodic witness against generate_collars(1).
#include "collars.hpp"
#include "hex_grid_svg.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
namespace h=jw::hexgrid;namespace g=jw::grid;
int main(int argc,char** argv) {
 try {
  if(argc!=3)throw std::invalid_argument("usage: verify INPUT OUTPUT_PREFIX");
  std::ifstream in(argv[1]);g::integer a,b,skew;in>>a>>b>>skew;
  if(!in || a<=0 || b<=0 || a*b>10000)throw std::invalid_argument("invalid period dimensions");
  h::hex_patch cell;g::integer q,r;std::string color;unsigned orientation;
  while(in>>q>>r>>color>>orientation) {
   if(q<0 || q>=a || r<0 || r>=b || !cell.insert({{q,r},g::parse_type(color),orientation}))throw std::invalid_argument("invalid cell");
  }
  if(cell.size()!=static_cast<size_t>(a*b))throw std::invalid_argument("incomplete period cell");
  auto wrap=[&](h::coordinate p) {
   auto quotient=p.r/b,rem=p.r%b;if(rem<0){rem+=b;--quotient;}
   auto u=g::subtract(p.q,g::multiply(quotient,skew))%a;if(u<0)u+=a;
   return h::coordinate{u,rem};
  };
  const auto catalog=h::generate_collars(1);std::set<h::collar> allowed(catalog.collars.begin(),catalog.collars.end());
  std::array<size_t,5> counts{};std::set<size_t> used;
  for(const auto& [at,tile]:cell.tiles) {
   h::collar local{{{},tile.type,tile.orientation}};++counts[unsigned(tile.type)];
   for(auto d:h::directions) {const auto& x=cell.at(wrap(at+d));local.push_back({d,x.type,x.orientation});}
   const auto key=h::canonicalize_collar(local);
   if(!allowed.contains(key))throw std::logic_error("wrapped collar absent from production catalog");
   used.insert(static_cast<size_t>(std::lower_bound(catalog.collars.begin(),catalog.collars.end(),key)-catalog.collars.begin()));
  }
  if(std::any_of(counts.begin(),counts.end(),[](size_t n){return !n;}))throw std::logic_error("missing color");
  h::to_svg(std::string(argv[2])+"-cell.svg",cell);
  h::hex_patch repeat;
  for(g::integer j=0;j<3;++j)for(g::integer i=0;i<3;++i)for(auto [at,tile]:cell.tiles) {
   tile.at=tile.at+h::coordinate{g::add(g::multiply(i,a),g::multiply(j,skew)),g::multiply(j,b)};
   repeat.insert(tile);
  }
  h::to_svg(std::string(argv[2])+"-repeated.svg",repeat);
  std::cout<<"Verified "<<cell.size()<<" wrapped collars against production catalog; period vectors ("<<a<<",0), ("<<skew<<','<<b<<").\n";
  for(size_t i=0;i<5;++i)std::cout<<g::name(g::types[i])<<": "<<counts[i]<<'\n';
  std::cout<<"Distinct radius-1 collar indices used:";for(auto i:used)std::cout<<' '<<i;std::cout<<'\n';
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
