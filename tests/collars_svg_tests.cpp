#include "collars.hpp"
#include "hex_grid_svg.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace h=jw::hexgrid;
void require(bool b,const char* s) {if(!b)throw std::runtime_error(s);}
template<class F> void fails(F f) {try {f();} catch(const std::exception&) {return;}throw std::runtime_error("expected SVG rejection");}
std::string read(const std::string& p) {std::ifstream in(p);return {std::istreambuf_iterator<char>(in),{}};}
int main() {
 try {
  std::filesystem::create_directories("collar-test-output");
  for(unsigned n:{0,1,6}) {
   const auto catalog=h::generate_collars(n);const auto before=catalog.collars;
   const std::string path="collar-test-output/collars-"+std::to_string(n)+".svg";
   h::collars_to_svg(path,catalog);const auto first=read(path);
   h::collars_to_svg(path,catalog);require(read(path)==first,"deterministic SVG");require(catalog.collars==before,"SVG mutates catalog");
  }
  h::collars_to_svg("collar-test-output/empty.svg",{});
  auto c=h::generate_collars(1);
  h::collars_to_svg("collar-test-output/single-column.svg",c,{1,4});
  fails([&]{h::collars_to_svg("collar-test-output/no-such-directory/fail.svg",c);});
  fails([&]{h::collars_to_svg("collar-test-output/invalid.svg",c,{0,12});});
  fails([&]{h::collars_to_svg("collar-test-output/invalid.svg",c,{2,std::numeric_limits<double>::infinity()});});
  c.collars[0].pop_back();fails([&]{h::collars_to_svg("collar-test-output/invalid.svg",c);});
  std::cout<<"SVG generation, determinism, const input and errors passed\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
