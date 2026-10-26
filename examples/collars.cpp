#include "collars.hpp"
#include "hex_grid_svg.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
int main(int argc,char** argv) {
    try {
        const auto number=argc>1?std::stoul(argv[1]):2ul;
        if(number>32) throw std::invalid_argument("radius must be <=32");
        const auto catalog=jw::hexgrid::generate_collars(static_cast<unsigned>(number));
        const std::string path=argc>2?argv[2]:"collars.svg";
        jw::hexgrid::collars_to_svg(path,catalog,{5,12});
        std::cout<<catalog.collars.size()<<" collars; depth "<<catalog.inflation_depth<<"; "<<path<<'\n';
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
