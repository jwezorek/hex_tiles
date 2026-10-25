#include "hex_grid_svg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
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
