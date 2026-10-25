#include "combinatorial_svg.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace jw::grid {
    namespace {
        struct point { double x,y; };
        point project(vertex v) {
            return {static_cast<double>(v.q)+static_cast<double>(v.r)/2.0,
                    std::sqrt(3.0)/2.0*static_cast<double>(v.r)};
        }
        struct bounds {
            double min_x=std::numeric_limits<double>::infinity(),min_y=min_x;
            double max_x=-min_x,max_y=-min_x;
            void include(point p) {
                if (!std::isfinite(p.x) || !std::isfinite(p.y)) throw std::overflow_error("nonfinite grid SVG coordinate");
                min_x=std::min(min_x,p.x); min_y=std::min(min_y,p.y);
                max_x=std::max(max_x,p.x); max_y=std::max(max_y,p.y);
            }
            double width() const { return max_x-min_x; }
            double height() const { return max_y-min_y; }
        };
        struct shape {
            tile_type type;
            node_id id,parent;
            std::size_t panel;
            std::vector<point> outline;
            std::vector<std::pair<point,point>> emphasized;
            point arrow_from{},arrow_to{};
        };

        std::vector<vertex> boundary(const occurrence& leaf) {
            // Every cell contributes a CCW directed boundary. Interior edges
            // occur twice with opposite directions and cancel exactly.
            std::map<edge,std::pair<vertex,vertex>> exposed;
            for (auto cell:leaf.cells) {
                const auto v=cell.vertices();
                for (std::size_t i=0;i<3;++i) {
                    const auto from=v[i],to=v[(i+1)%3];
                    auto [it,inserted]=exposed.emplace(edge{from,to},std::pair{from,to});
                    if (!inserted) {
                        if (it->second!=std::pair{to,from}) throw std::logic_error("overlapping leaf cells in SVG boundary");
                        exposed.erase(it);
                    }
                }
            }
            std::map<vertex,vertex> next;
            for (const auto& [edge,segment]:exposed) {
                (void)edge;
                if (!next.emplace(segment).second) throw std::logic_error("branched leaf boundary");
            }
            if (next.empty()) throw std::invalid_argument("cannot render an empty leaf");
            std::vector<vertex> result;
            const vertex start=next.begin()->first;
            vertex current=start;
            do {
                result.push_back(current);
                auto it=next.find(current);
                if (it==next.end()) throw std::logic_error("open leaf boundary");
                current=it->second;
                next.erase(it);
            } while (current!=start);
            if (!next.empty()) throw std::logic_error("disconnected leaf boundary");
            return result;
        }
        shape draw_leaf(const occurrence& leaf,node_id id,std::size_t panel,bool arrows) {
            shape result{leaf.type,id,leaf.parent,panel,{},{},{},{}};
            const auto outline=boundary(leaf);
            const auto inverse=leaf.transform.inverse();
            const bool triangle=leaf.type==tile_type::yellow || leaf.type==tile_type::red;
            bool have_base=false;
            vertex left{},right{};
            for (std::size_t i=0;i<outline.size();++i) {
                const auto from=outline[i],to=outline[(i+1)%outline.size()];
                const auto u=inverse(from),v=inverse(to);
                result.outline.push_back(project(from));
                // Decoration on the exposed CELL edges, matching external_edges
                // in hex.cpp. This is distinct from the whole patch's perimeter.
                const bool bold=triangle ? (add(u.q,u.r)==1 && add(v.q,v.r)==1)
                                         : !(u.r==0 && v.r==0);
                if (bold) result.emphasized.emplace_back(project(from),project(to));
                if (u.r==0) {
                    if (!have_base) { left=right=u; have_base=true; }
                    if (u.q<left.q) left=u;
                    if (u.q>right.q) right=u;
                }
            }
            if (arrows) {
                if (!have_base || left==right) throw std::logic_error("missing lattice arrow endpoints");
                const auto a=project(leaf.transform(left)),b=project(leaf.transform(right));
                const double dx=b.x-a.x,dy=b.y-a.y;
                // Reflecting an inset arrow also reflects its normal. Offsetting
                // every arrow to the world-left would be wrong for mirrored leaves.
                const double sign=static_cast<double>(leaf.transform.determinant());
                const point normal{-0.06*dy*sign,0.06*dx*sign};
                result.arrow_from={a.x+0.12*dx+normal.x,a.y+0.12*dy+normal.y};
                result.arrow_to={b.x-0.12*dx+normal.x,b.y-0.12*dy+normal.y};
            }
            return result;
        }
        std::vector<shape> draw_patch(const patch& p,std::size_t panel,bool arrows) {
            std::vector<shape> result;
            result.reserve(p.leaves.size());
            for (auto id:p.leaves) result.push_back(draw_leaf(p.nodes.at(id),id,panel,arrows));
            return result;
        }
        bounds measure(const std::vector<shape>& shapes) {
            bounds result;
            for (const auto& s:shapes) for (auto v:s.outline) result.include(v);
            return result;
        }
        void move(shape& s,double dx,double dy) {
            const auto shift=[&](point& p) { p.x+=dx; p.y+=dy; };
            for (auto& p:s.outline) shift(p);
            for (auto& [a,b]:s.emphasized) { shift(a); shift(b); }
            shift(s.arrow_from); shift(s.arrow_to);
        }
        void write(const std::string& path,const std::vector<shape>& shapes,bool arrows,double extent) {
            if (shapes.empty()) throw std::invalid_argument("cannot render an empty grid patch");
            const auto b=measure(shapes);
            const double span=std::max(b.width(),b.height());
            if (!(span>0) || !std::isfinite(span)) throw std::overflow_error("invalid grid SVG bounds");
            const double scale=extent/span,margin=20.0;
            const double width=b.width()*scale+2*margin,height=b.height()*scale+2*margin;
            const auto svg=[&](point p) {
                point result{(p.x-b.min_x)*scale+margin,(b.max_y-p.y)*scale+margin};
                if (!std::isfinite(result.x) || !std::isfinite(result.y)) throw std::overflow_error("nonfinite grid SVG projection");
                return result;
            };
            std::ofstream out(path);
            if (!out) throw std::runtime_error("could not open grid SVG: "+path);
            out<<std::setprecision(17)<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
               <<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<width<<"\" height=\""<<height
               <<"\" viewBox=\"0 0 "<<width<<' '<<height<<"\">\n";
            if (arrows) out<<"<defs><marker id=\"arrowhead\" markerWidth=\"4\" markerHeight=\"3\" refX=\"3.8\" refY=\"1.5\" orient=\"auto\" markerUnits=\"strokeWidth\"><polygon points=\"0 0, 4 1.5, 0 3\" fill=\"#000000\" /></marker></defs>\n";
            out<<"<g id=\"subtiles\">\n";
            for (const auto& s:shapes) {
                out<<"<g class=\"subtile\" data-type=\""<<name(s.type)<<"\" data-panel=\""<<s.panel
                   <<"\" data-node=\""<<s.id<<'"';
                if (s.parent!=no_parent) out<<" data-parent=\""<<s.parent<<'"';
                out<<">\n<polygon points=\"";
                for (auto p:s.outline) { const auto v=svg(p); out<<v.x<<','<<v.y<<' '; }
                out<<"\" fill=\""<<color(s.type)<<"\" stroke=\"none\" />\n";
                out<<"<path d=\"";
                for (const auto& [a,b]:s.emphasized) {
                    const auto u=svg(a),v=svg(b);
                    out<<"M "<<u.x<<' '<<u.y<<" L "<<v.x<<' '<<v.y<<' ';
                }
                out<<"\" fill=\"none\" stroke=\"#000000\" stroke-width=\"1\" stroke-linecap=\"round\" stroke-linejoin=\"round\" />\n";
                if (arrows) {
                    const auto a=svg(s.arrow_from),b=svg(s.arrow_to);
                    out<<"<line x1=\""<<a.x<<"\" y1=\""<<a.y<<"\" x2=\""<<b.x<<"\" y2=\""<<b.y
                       <<"\" stroke=\"#000000\" stroke-width=\"1.5\" stroke-linecap=\"round\" marker-end=\"url(#arrowhead)\" />\n";
                }
                out<<"</g>\n";
            }
            out<<"</g>\n</svg>\n";
            out.close();
            if (!out) throw std::runtime_error("writing grid SVG failed: "+path);
        }
    }
    void to_svg(const std::string& path,const patch& value,bool show_arrows) {
        write(path,draw_patch(value,0,show_arrows),show_arrows,1200.0);
    }
    void rules_svg(const std::string& path,bool show_arrows) {
        std::vector<shape> diagram;
        double cursor=0;
        std::size_t panel=0;
        for (auto type:types) {
            auto leaf=draw_patch(make_supertile(type,0),panel++,show_arrows);
            auto children=draw_patch(make_supertile(type,1),panel++,show_arrows);
            const auto lb=measure(leaf),cb=measure(children);
            const double row_height=std::max(lb.height(),cb.height());
            const auto append=[&](std::vector<shape>& source,const bounds& b,double left) {
                const double dx=left-b.min_x,dy=-cursor-(row_height-b.height())/2-b.max_y;
                for (auto& s:source) { move(s,dx,dy); diagram.push_back(std::move(s)); }
            };
            append(leaf,lb,0);
            append(children,cb,lb.width()+1.2);
            cursor+=row_height+0.9;
        }
        write(path,diagram,show_arrows,1300.0);
    }
}
