#include "hex_grid_svg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>

namespace jw::hexgrid {
    namespace {
        struct point { double x=0,y=0; };
        point project(grid::vertex p) {
            return {static_cast<double>(p.q)+0.5*static_cast<double>(p.r),
                    std::sqrt(3.0)/2.0*static_cast<double>(p.r)};
        }
        struct bounds {
            double left=std::numeric_limits<double>::infinity();
            double right=-std::numeric_limits<double>::infinity();
            double bottom=left,top=right;
            void include(point p) {
                if (!std::isfinite(p.x) || !std::isfinite(p.y))
                    throw std::overflow_error("nonfinite hex SVG coordinate");
                left=std::min(left,p.x); right=std::max(right,p.x);
                bottom=std::min(bottom,p.y); top=std::max(top,p.y);
            }
        };
        constexpr std::array<grid::vertex,6> corners{{
            {1,0},{0,1},{-1,1},{-1,0},{0,-1},{1,-1}
        }};
    }
    void to_svg(const std::string& path,const hex_patch& patch) {
        bounds b;
        for (const auto& [at,tile]:patch.tiles) {
            const auto center=patch.placement.center(at);
            if (tile.at!=at || tile.orientation>=6)
                throw std::logic_error("invalid hex tile coordinate or orientation");
            for (auto corner:corners) b.include(project(center+corner));
        }
        // An empty patch (e.g. depth zero) has no complete hexagons.
        if (patch.empty()) { b.left=b.bottom=0; b.right=b.top=1; }
        const double span=std::max(b.right-b.left,b.top-b.bottom);
        if (!(span>0) || !std::isfinite(span)) throw std::overflow_error("invalid hex SVG bounds");
        const double scale=1200.0/span,margin=20.0;
        const double width=(b.right-b.left)*scale+2*margin;
        const double height=(b.top-b.bottom)*scale+2*margin;
        const auto svg=[&](point p)->point {
            point v{(p.x-b.left)*scale+margin,(b.top-p.y)*scale+margin};
            if (!std::isfinite(v.x) || !std::isfinite(v.y))
                throw std::overflow_error("nonfinite projected hex SVG coordinate");
            return v;
        };
        std::ofstream out(path);
        if (!out) throw std::runtime_error("could not open hex SVG: "+path);
        out<<std::setprecision(17)
           <<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           <<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<width
           <<"\" height=\""<<height<<"\" viewBox=\"0 0 "<<width<<' '<<height<<"\">\n"
           <<"<defs><marker id=\"hex-arrowhead\" markerWidth=\"4\" markerHeight=\"3\" "
           <<"refX=\"3.8\" refY=\"1.5\" orient=\"auto\" markerUnits=\"strokeWidth\">"
           <<"<polygon points=\"0 0, 4 1.5, 0 3\" fill=\"#151515\" /></marker></defs>\n"
           <<"<g id=\"hexagons\">\n";
        for (const auto& [at,tile]:patch.tiles) {
            const auto center=patch.placement.center(at);
            const auto middle=project(center);
            const auto start=svg(middle);
            const double angle=std::numbers::pi/3.0*static_cast<double>(tile.orientation);
            const point tip=svg({middle.x+0.54*std::cos(angle),
                                 middle.y+0.54*std::sin(angle)});
            out<<"<g class=\"hex-tile\" data-type=\""<<grid::name(tile.type)
               <<"\" data-q=\""<<at.q<<"\" data-r=\""<<at.r
               <<"\" data-orientation=\""<<tile.orientation<<"\">\n"
               <<"<polygon points=\"";
            for (auto corner:corners) {
                const auto p=svg(project(center+corner));
                out<<p.x<<','<<p.y<<' ';
            }
            out<<"\" fill=\""<<grid::color(tile.type)
               <<"\" stroke=\"#202020\" stroke-width=\"1\" stroke-linejoin=\"round\" />\n"
               <<"<line x1=\""<<start.x<<"\" y1=\""<<start.y
               <<"\" x2=\""<<tip.x<<"\" y2=\""<<tip.y
               <<"\" stroke=\"#151515\" stroke-width=\"1.3\" stroke-linecap=\"round\" "
               <<"marker-end=\"url(#hex-arrowhead)\" />\n</g>\n";
        }
        out<<"</g>\n</svg>\n";
        out.close();
        if (!out) throw std::runtime_error("writing hex SVG failed: "+path);
    }
}

namespace jw::hexgrid {
    void collars_to_svg(const std::string& path,const collar_catalog& catalog,const collar_svg_options& options) {
        if(options.columns==0 || options.columns>1000 || !std::isfinite(options.hex_size) ||
           options.hex_size<4 || options.hex_size>1000 || catalog.radius>32)
            throw std::invalid_argument("invalid collar SVG radius, columns or hex size");
        const auto n=grid::integer(catalog.radius);
        const auto expected=grid::add(1,grid::multiply(3,grid::multiply(n,grid::add(n,1))));
        if(catalog.collars.size()>1'000'000/static_cast<size_t>(expected))
            throw std::length_error("collar SVG exceeds 1,000,000-hex rendering budget");
        bounds b;
        for(const auto& collar:catalog.collars) {
            if(collar.size()!=static_cast<size_t>(expected)) throw std::invalid_argument("incomplete collar in SVG catalog");
            std::set<coordinate> positions;
            for(auto tile:collar) {
                const bool small=tile.type==grid::tile_type::yellow || tile.type==grid::tile_type::red;
                if(static_cast<unsigned>(tile.type)>=grid::types.size() || tile.orientation>=(small?2u:6u) ||
                   tile.at.q < -n || tile.at.q>n || tile.at.r < -n || tile.at.r>n ||
                   grid::add(tile.at.q,tile.at.r)<-n || grid::add(tile.at.q,tile.at.r)>n || !positions.insert(tile.at).second)
                    throw std::invalid_argument("invalid collar hex in SVG catalog");
                const auto center=lattice{}.center(tile.at);
                for(auto corner:corners) b.include(project(center+corner));
            }
        }
        if(catalog.collars.empty()) {b.left=b.bottom=0;b.right=b.top=1;}
        const double size=options.hex_size,margin=16,gap=24,header=44,label=24;
        const auto cols=std::max(size_t(1),std::min(options.columns,catalog.collars.size()));
        const auto rows=catalog.collars.empty()?size_t(1):1+(catalog.collars.size()-1)/cols;
        const double panel_w=(b.right-b.left)*size+2*margin;
        const double panel_h=(b.top-b.bottom)*size+2*margin+label;
        const double width=std::max(300.0,cols*panel_w+(cols+1)*gap);
        const double height=header+rows*panel_h+(rows+1)*gap;
        if(!std::isfinite(width) || !std::isfinite(height)) throw std::overflow_error("invalid collar SVG bounds");
        std::ofstream out(path);
        if(!out) throw std::runtime_error("could not open collar SVG: "+path);
        out<<std::setprecision(17)<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           <<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<width<<"\" height=\""<<height
           <<"\" viewBox=\"0 0 "<<width<<' '<<height<<"\" data-radius=\""<<catalog.radius
           <<"\" data-count=\""<<catalog.collars.size()<<"\" data-inflation-depth=\""<<catalog.inflation_depth
           <<"\" data-depth-basis=\""<<(catalog.depth_basis==collar_depth_basis::certified_catalog?"certified-catalog":
                                      catalog.depth_basis==collar_depth_basis::coverage_certificate?"coverage-certificate":"single-colors")<<"\">\n"
           <<"<title>Radius "<<catalog.radius<<" collars ("<<catalog.collars.size()<<")</title>\n"
           <<"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n"
           <<"<defs><marker id=\"collar-arrowhead\" markerWidth=\""<<size*.20<<"\" markerHeight=\""<<size*.18
           <<"\" refX=\"0.9\" refY=\"0.5\" viewBox=\"0 0 1 1\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">"
           <<"<path d=\"M0,0 L1,0.5 L0,1 Z\" fill=\"#151515\"/></marker></defs>\n"
           <<"<text x=\"24\" y=\"28\" font-family=\"sans-serif\" font-size=\"20\">Radius "<<catalog.radius
           <<" · "<<catalog.collars.size()<<" collars</text>\n";
        if(catalog.collars.empty()) out<<"<text x=\"24\" y=\"70\" font-family=\"sans-serif\" font-size=\"14\">Empty catalog</text>\n";
        for(size_t i=0;i<catalog.collars.size();++i) {
            const double left=gap+(i%cols)*(panel_w+gap),top=header+gap+(i/cols)*(panel_h+gap);
            auto svg=[&](point p) {return point{left+margin+(p.x-b.left)*size,top+label+margin+(b.top-p.y)*size};};
            out<<"<g class=\"collar-panel\" id=\"collar-"<<i<<"\" data-collar-index=\""<<i<<"\">\n"
               <<"<text x=\""<<left+margin<<"\" y=\""<<top+16<<"\" font-family=\"sans-serif\" font-size=\"14\">#"<<i<<"</text>\n";
            for(auto tile:catalog.collars[i]) {
                const auto center=lattice{}.center(tile.at);const auto p=project(center);const auto start=svg(p);
                const auto angle=std::numbers::pi/3.0*tile.orientation;
                const auto tip=svg({p.x+.54*std::cos(angle),p.y+.54*std::sin(angle)});
                const bool central=tile.at==coordinate{};
                out<<"<g class=\"hex-tile\" data-collar-index=\""<<i<<"\" data-q=\""<<tile.at.q
                   <<"\" data-r=\""<<tile.at.r<<"\" data-color=\""<<grid::name(tile.type)<<"\" data-orientation=\""<<tile.orientation
                   <<"\" data-central=\""<<(central?"true":"false")<<"\">\n<polygon points=\"";
                for(auto corner:corners) {const auto v=svg(project(center+corner));out<<v.x<<','<<v.y<<' ';}
                out<<"\" fill=\""<<grid::color(tile.type)<<"\" stroke=\""<<(central?"#000000":"#555555")
                   <<"\" stroke-width=\""<<(central?3:0.7)<<"\" stroke-linejoin=\"round\"/>\n"
                   <<"<line class=\"orientation-arrow\" x1=\""<<start.x<<"\" y1=\""<<start.y<<"\" x2=\""<<tip.x
                   <<"\" y2=\""<<tip.y<<"\" stroke=\"#151515\" stroke-width=\""<<size*.065
                   <<"\" marker-end=\"url(#collar-arrowhead)\"/>\n</g>\n";
            }
            out<<"</g>\n";
        }
        out<<"</svg>\n";out.close();
        if(!out) throw std::runtime_error("writing collar SVG failed: "+path);
    }
}
