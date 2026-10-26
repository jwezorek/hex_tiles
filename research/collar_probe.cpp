// Research probe, deliberately separate from the production converter/enumerator.
// Build instructions and interpretation are in collar-depth.md.
#include "finite_patch.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace g=jw::grid;
namespace h=jw::hexgrid;
using P=g::vertex;
using C=h::coordinate;
using Key=std::vector<std::tuple<long long,long long,int,unsigned>>;
static int mod(int a,int b) {return (a%b+b)%b;}
static C rot(C a) {return {-a.r,a.q+a.r};}
static int dist(C a,C b={}) {auto d=a-b;return int(std::max({std::abs(d.q),std::abs(d.r),std::abs(d.q+d.r)}));}
static std::vector<C> ball(int n) {
    std::vector<C> b;
    for(int q=-n;q<=n;++q) for(int r=-n;r<=n;++r)
        if(dist({q,r})<=n) b.push_back({q,r});
    return b;
}
static Key key(const h::hex_patch& p) {
    std::optional<Key> best;
    for(unsigned t=0;t<6;++t) {
        Key k;
        for(auto [pos,x]:p.tiles) {
            auto at=pos;
            for(unsigned j=0;j<t;++j) at=rot(at);
            const unsigned m=(x.type==g::tile_type::yellow || x.type==g::tile_type::red)?2:6;
            k.emplace_back(at.q,at.r,int(x.type),(x.orientation+t)%m);
        }
        std::sort(k.begin(),k.end());
        auto [q,r,c,o]=k.front();
        for(auto& [x,y,z,w]:k) {x-=q;y-=r;}
        if(!best || k<*best) best=k;
    }
    return *best;
}
static std::vector<P> corners(g::tile_type t,g::lengths d) {
    bool tiny=t==g::tile_type::yellow || t==g::tile_type::red;
    auto z=(t==g::tile_type::purple || t==g::tile_type::red)?d.S:d.middle();
    if(tiny) return {{0,0},{z,0},{0,z}};
    return {{0,0},{d.L,0},{d.L-z,z},{0,z}};
}
static std::map<P,P> deform(const g::patch& src,g::lengths d) {
    struct Item {size_t id;P local;};
    std::map<P,std::vector<Item>> junction;
    for(auto id:src.leaves) {
        auto& x=src.nodes[id];auto old=corners(x.type,{}),nw=corners(x.type,d);
        auto lin=x.transform;lin.translation={};
        for(size_t i=0;i<old.size();++i) junction[x.transform(old[i])].push_back({id,lin(nw[i])});
    }
    std::vector<std::vector<std::pair<size_t,P>>> adj(src.nodes.size());
    for(auto& [v,a]:junction) for(size_t i=1;i<a.size();++i) {
        auto delta=a[0].local-a[i].local;
        adj[a[0].id].push_back({a[i].id,delta});adj[a[i].id].push_back({a[0].id,P{}-delta});
    }
    std::vector<std::optional<P>> offset(src.nodes.size());
    size_t root=src.leaves.front();offset[root]=src.nodes[root].transform.translation;
    std::queue<size_t> q;q.push(root);
    while(!q.empty()) {auto at=q.front();q.pop();for(auto [to,diff]:adj[at]) {
        auto v=*offset[at]+diff;
        if(!offset[to]) {offset[to]=v;q.push(to);} else if(*offset[to]!=v) throw std::logic_error("corner inconsistency");
    }}
    std::map<P,P> out;
    for(auto& [v,a]:junction) {if(!offset[a[0].id]) throw std::logic_error("disconnected mesh");out[v]=*offset[a[0].id]+a[0].local;}
    return out;
}
static constexpr std::array<P,6> unit{{{1,0},{0,1},{-1,1},{-1,0},{0,-1},{1,-1}}};
static std::array<g::cell,6> sectors(P p) {
    return {{{p+P{-1,-1},true},{p+P{0,-1},false},{p+P{0,-1},true},
             {p,false},{p+P{-1,0},true},{p+P{-1,0},false}}};
}
// Assumed improved reconstruction, implemented only in this research probe.
// Any observed sector uniquely selects one of the repository's motifs.
static h::hex_patch partial(const g::patch& p,h::lattice lattice) {
    h::hex_patch out;out.placement=lattice;
    for(auto [cell,id]:p.occupancy) for(auto v:cell.vertices()) {
        C at;
        try {at=lattice.locate(v);} catch(const std::invalid_argument&) {continue;}
        auto ring=sectors(v);int s=int(std::find(ring.begin(),ring.end(),cell)-ring.begin());
        auto& leaf=p.nodes[id];P arrow{leaf.transform.a,leaf.transform.c};
        int dir=int(std::find(unit.begin(),unit.end(),arrow)-unit.begin());
        bool mirror=leaf.transform.determinant()<0;
        bool tiny=leaf.type==g::tile_type::yellow || leaf.type==g::tile_type::red;
        std::vector<unsigned> matches;
        for(int o=0;o<(tiny?2:6);++o) {
            int j=mod(s-o,6);
            int expected=tiny?(j<2?4:j<4?0:2):0;
            bool reflected=tiny?(j%2==0):(j<3);
            if(mod(expected+o,6)==dir && reflected==mirror) matches.push_back(unsigned(o));
        }
        if(matches.size()!=1) throw std::logic_error("one sector failed unique motif decoding");
        h::hex_tile x{at,leaf.type,matches[0]};
        auto it=out.tiles.find(at);
        if(it!=out.tiles.end() && (it->second.type!=x.type || it->second.orientation!=x.orientation))
            throw std::logic_error("inconsistent partial motif");
        out.insert(x);
    }
    return out;
}
static std::map<P,std::vector<C>> vertices(const h::hex_patch& p) {
    std::map<P,std::vector<C>> out;
    for(auto& [at,x]:p.tiles) for(auto u:unit) out[p.placement.center(at)+u].push_back(at);
    return out;
}
static std::map<Key,h::hex_patch> sample_vertices(unsigned depth) {
    std::map<Key,h::hex_patch> out;
    for(auto type:g::types) {
        auto hex=h::convert(g::make_supertile(type,depth));
        for(auto& [v,cs]:vertices(hex)) if(cs.size()==3) {
            h::hex_patch star;star.placement=hex.placement;
            for(auto at:cs) star.insert(hex.at(at));
            auto k=key(star);
            if(!out.contains(k)) {
                // Normalize translations so every sample starts close to zero.
                h::hex_patch norm;auto first=star.tiles.begin()->first;
                for(auto& [at,x]:star.tiles) {auto y=x;y.at=at-first;norm.insert(y);}
                out.emplace(k,std::move(norm));
            }
        }
    }
    return out;
}

static std::map<Key,h::hex_patch> sample_balls(unsigned depth,int radius) {
    std::map<Key,h::hex_patch> out;auto b=ball(radius);
    for(auto type:g::types) {
        auto hex=h::convert(g::make_supertile(type,depth));
        for(auto& [center,tile]:hex.tiles) {
            h::hex_patch p;
            for(auto d:b) {auto x=hex.find(center+d);if(!x) break;auto y=*x;y.at=d;p.insert(y);}
            if(p.size()==b.size()) out.try_emplace(key(p),std::move(p));
        }
    }
    return out;
}
static std::map<Key,h::hex_patch> compatible_balls1(const std::map<Key,h::hex_patch>& allowed) {
    std::vector<h::hex_tile> labels;
    for(auto t:g::types) for(unsigned o=0;o<((t==g::tile_type::yellow || t==g::tile_type::red)?2u:6u);++o)
        labels.push_back({{},t,o});
    std::map<Key,h::hex_patch> out;
    for(auto center:labels) {
        std::array<h::hex_tile,6> ring;
        auto ok=[&](int edge) {
            h::hex_patch p;p.insert(center);p.insert(ring[edge]);p.insert(ring[(edge+1)%6]);
            return allowed.contains(key(p));
        };
        auto visit=[&](auto&& self,int j)->void {
            if(j==6) {
                if(!ok(5)) return;
                h::hex_patch p;p.insert(center);for(auto x:ring) p.insert(x);out.try_emplace(key(p),p);return;
            }
            for(auto x:labels) {x.at=h::directions[j];ring[j]=x;if(j && !ok(j-1)) continue;self(self,j+1);}
        };
        visit(visit,0);
    }
    return out;
}
static std::map<Key,h::hex_patch> compatible_balls2(const std::map<Key,h::hex_patch>& allowed1) {
    std::map<Key,h::hex_patch> variants;
    for(auto& [sig,p]:allowed1) for(unsigned t=0;t<6;++t) {
        h::hex_patch r;
        for(auto [at,tile]:p.tiles) {
            auto x=tile;for(unsigned j=0;j<t;++j) x.at=rot(x.at);
            unsigned m=(x.type==g::tile_type::yellow || x.type==g::tile_type::red)?2:6;
            x.orientation=(x.orientation+t)%m;r.insert(x);
        }
        Key raw;for(auto [at,x]:r.tiles) raw.emplace_back(at.q,at.r,int(x.type),x.orientation);
        variants.try_emplace(raw,r);
    }
    std::map<Key,h::hex_patch> out;
    auto visit=[&](auto&& self,h::hex_patch p,int j)->void {
        if(j==6) {if(p.size()!=19) throw std::logic_error("radius2 join has wrong size");out.try_emplace(key(p),p);return;}
        C center=h::directions[j];
        for(auto& [sig,templ]:variants) {
            bool ok=true;
            for(auto [at,x]:templ.tiles) {auto y=p.find(center+at);if(y && (x.type!=y->type || x.orientation!=y->orientation)) {ok=false;break;}}
            if(!ok) continue;
            auto next=p;for(auto [at,x]:templ.tiles) {x.at=center+at;next.insert(x);}self(self,std::move(next),j+1);
        }
    };
    for(auto& [sig,p]:variants) visit(visit,p,0);
    return out;
}
using Polygon=std::array<P,6>;
static bool inside3(P triple,const Polygon& polygon) {
    for(size_t i=0;i<6;++i) {
        P a=polygon[i],d=polygon[(i+1)%6]-a;
        auto x=triple-P{3*a.q,3*a.r};
        if(d.q*x.r-d.r*x.q<0) return false;
    }
    return true;
}
static std::set<C> footprint(const Polygon& polygon,h::lattice lattice) {
    long long qmin=polygon[0].q,qmax=qmin,rmin=polygon[0].r,rmax=rmin;
    for(auto p:polygon) {qmin=std::min(qmin,p.q);qmax=std::max(qmax,p.q);rmin=std::min(rmin,p.r);rmax=std::max(rmax,p.r);}
    std::set<C> out;
    for(auto q=qmin-1;q<=qmax+1;++q) for(auto r=rmin-1;r<=rmax+1;++r) {
        P p{q,r};C at;
        try {at=lattice.locate(p);} catch(const std::invalid_argument&) {continue;}
        for(auto sector:sectors(p)) {
            P sum{};for(auto v:sector.vertices()) sum=sum+v;
            if(inside3(sum,polygon)) {out.insert(at);break;}
        }
    }
    return out;
}
static h::lattice direct_phase(const g::patch& src,const std::map<P,P>& moved,unsigned k) {
    if(!k) return {};
    // A complete base prototile already fixes its hex center: canonical (0,0)
    // for a tiny triangle, (1,0) for a trapezoid. Follow one descendant only.
    const auto& first=src.nodes[src.leaves.front()];auto frame=first.transform;
    frame.translation=moved.at(first.transform(P{0,0}));auto type=first.type;
    std::vector<g::lengths> levels;g::lengths d;
    for(unsigned j=0;j<k;++j) {levels.push_back(d);d=g::next_lengths(d);}
    for(unsigned j=k;j>0;--j) {
        auto child=g::rules(levels[j-1])[unsigned(type)].front();frame=frame*child.transform;type=child.type;
    }
    bool tiny=type==g::tile_type::yellow || type==g::tile_type::red;
    auto center=frame(tiny?P{}:P{1,0});
    return h::lattice{{mod(int(center.q%3-center.r%3),3),0}};
}
struct GeometryResult {int radius=999999;C target{};h::lattice phase;Polygon central;std::vector<std::pair<P,C>> missing;};
static GeometryResult geometry(const h::hex_patch& source,unsigned k) {
    auto src=h::to_triangles(source);g::lengths dims;for(unsigned i=0;i<k;++i) dims=g::next_lengths(dims);
    auto moved=deform(src,dims);auto phase=k?direct_phase(src,moved,k):source.placement;
    std::map<C,Polygon> polygons;std::map<C,std::set<C>> images;
    for(auto& [at,tile]:source.tiles) {
        Polygon poly;for(size_t i=0;i<6;++i) poly[i]=moved.at(source.placement.center(at)+unit[i]);
        polygons[at]=poly;images[at]=footprint(poly,phase);
    }
    auto allvertices=vertices(source);
    std::map<C,int> best;
    std::map<P,std::set<C>> stars;
    for(auto& [v,cs]:allvertices) if(cs.size()==3) {
        auto& image=stars[v];for(auto a:cs) image.insert(images[a].begin(),images[a].end());
        std::map<C,int> distance;std::queue<C> todo;
        for(auto x:image) {
            auto ns=h::neighbors(x);
            if(std::any_of(ns.begin(),ns.end(),[&](C y){return !image.contains(y);})) {distance[x]=0;todo.push(x);}
        }
        while(!todo.empty()) {auto x=todo.front();todo.pop();for(auto y:h::neighbors(x)) {
            if(image.contains(y) && !distance.contains(y)) {distance[y]=distance[x]+1;todo.push(y);}
        }}
        for(auto [x,d]:distance) {
            auto found=best.find(x);if(found==best.end()) best.emplace(x,d);else found->second=std::max(found->second,d);
        }
    }
    GeometryResult result;result.phase=phase;result.central=polygons.at({});
    for(auto x:images.at({})) {
        P p=phase.center(x);
        if(!inside3({3*p.q,3*p.r},polygons.at({}))) continue;
        int r=best.contains(x)?best[x]:-1;
        if(r<result.radius) {result.radius=r;result.target=x;}
    }
    for(auto& [v,image]:stars) if(image.contains(result.target)) {
        for(auto d:ball(result.radius+1)) if(!image.contains(result.target+d)) {
            result.missing.push_back({v,result.target+d});break;
        }
    }
    return result;
}

int main(int argc,char** argv) {
 try {
    std::string mode=argc>1?argv[1]:"vertices";
    unsigned depth=argc>2?unsigned(std::stoul(argv[2])):4;
    if(mode=="extendibility") {
        auto p=g::make_supertile(g::tile_type::green,depth);g::lengths dims;
        for(unsigned j=0;j<depth;++j) dims=g::next_lengths(dims);
        std::set<g::tile_type> found;
        for(auto id:p.leaves) {
            const auto& n=p.nodes[id];if(n.transform.determinant()!=1 || found.contains(n.type)) continue;
            bool strict=true;
            for(auto c:n.cells) for(auto v:c.vertices())
                if(v.q<=0 || v.r<=0 || v.q+v.r>=dims.L || v.r>=dims.middle()) strict=false;
            if(strict) {found.insert(n.type);std::cout<<"interior_seed,"<<g::name(n.type)<<",leaf,"<<id<<",frame,"
                <<n.transform.a<<','<<n.transform.b<<','<<n.transform.c<<','<<n.transform.d<<','
                <<n.transform.translation.q<<','<<n.transform.translation.r<<'\n';}
        }
        if(found.size()!=5) throw std::logic_error("not all seed colors have strict interior occurrences");
        std::cout<<"all_five_seeds_extend_via_green_depth,"<<depth<<'\n';return 0;
    }
    if(mode=="vertices") {
        auto v=sample_vertices(depth);std::cout<<"depth,"<<depth<<",vertex_types,"<<v.size()<<'\n';
        return 0;
    }
    if(mode=="contexts") {
        std::cout<<"# type IDs: 0 yellow, 1 blue, 2 green, 3 purple, 4 red\n# entries q,r,type,orientation; radius2 centers are (0,0)\n";
        for(auto& [sig,p]:sample_vertices(depth)) {
            std::cout<<"vertex";for(auto [at,x]:p.tiles) std::cout<<';'<<at.q<<','<<at.r<<','<<int(x.type)<<','<<x.orientation;std::cout<<'\n';
        }
        for(auto& [sig,p]:sample_balls(depth,2)) {
            std::cout<<"radius2";for(auto [at,x]:p.tiles) std::cout<<';'<<at.q<<','<<at.r<<','<<int(x.type)<<','<<x.orientation;std::cout<<'\n';
        }
        return 0;
    }
    if(mode=="selfcheck") {
        size_t sectors_checked=0;
        for(auto t:g::types) for(unsigned o=0;o<((t==g::tile_type::yellow || t==g::tile_type::red)?2u:6u);++o) {
            h::hex_patch one;one.insert({{},t,o});auto full=h::to_triangles(one);
            for(auto [cell,id]:full.occupancy) {
                auto fragment=full;fragment.occupancy.clear();fragment.occupancy.emplace(cell,id);
                auto decoded=partial(fragment,{});
                if(decoded.size()!=1 || decoded.at({}).type!=t || decoded.at({}).orientation!=o)
                    throw std::logic_error("single-sector decoding mismatch");
                ++sectors_checked;
            }
        }
        size_t geometry_checked=0;
        for(auto& [sig,star]:sample_vertices(depth)) for(unsigned k=1;k<=3;++k) {
            auto src=h::to_triangles(star);g::lengths dims;for(unsigned j=0;j<k;++j) dims=g::next_lengths(dims);
            auto moved=deform(src,dims);auto phase=direct_phase(src,moved,k);
            g::patch direct;
            for(auto id:src.leaves) {
                const auto& n=src.nodes[id];auto parent=n.transform;parent.translation=moved.at(n.transform(P{}));
                auto super=g::make_supertile(n.type,k);
                for(auto childid:super.leaves) {
                    auto x=super.nodes[childid];x.parent=g::no_parent;x.transform=parent*x.transform;
                    for(auto& c:x.cells) c=parent(c);
                    auto newid=direct.nodes.size();direct.nodes.push_back(x);direct.leaves.push_back(newid);
                    for(auto c:x.cells) if(!direct.occupancy.emplace(c,newid).second) throw std::logic_error("direct overlap");
                }
            }
            auto actual=h::substitute_triangles(src,k);
            if(actual.leaves.size()!=direct.leaves.size()) throw std::logic_error("direct/iterated leaf count");
            auto offset=direct.nodes.front().transform.translation-actual.nodes.front().transform.translation;
            for(size_t j=0;j<direct.nodes.size();++j) {
                auto frame=actual.nodes[j].transform;frame.translation=frame.translation+offset;
                if(frame!=direct.nodes[j].transform || actual.nodes[j].type!=direct.nodes[j].type)
                    throw std::logic_error("direct/iterated rule placement");
            }
            std::set<C> expected;
            for(auto& [at,x]:star.tiles) {
                Polygon poly;for(size_t j=0;j<6;++j) poly[j]=moved.at(star.placement.center(at)+unit[j]);
                auto image=footprint(poly,phase);expected.insert(image.begin(),image.end());
            }
            auto decoded=partial(direct,phase);
            if(decoded.size()!=expected.size()) throw std::logic_error("polygon/triangle footprint count");
            for(auto [at,x]:decoded.tiles) if(!expected.contains(at)) throw std::logic_error("polygon/triangle footprint coordinate");
            ++geometry_checked;
        }
        std::cout<<"single_sector_checks,"<<sectors_checked<<",direct_vs_iterated_geometry_checks,"<<geometry_checked<<'\n';return 0;
    }
    if(mode=="admissible") {
        std::vector<h::hex_tile> labels;
        for(auto t:g::types) for(unsigned o=0;o<((t==g::tile_type::yellow || t==g::tile_type::red)?2u:6u);++o)
            labels.push_back({{},t,o});
        std::map<Key,h::hex_patch> accepted;
        size_t tried=0,geometric=0;auto known=sample_vertices(depth);
        for(auto a:labels) for(auto b:labels) for(auto c:labels) {
            a.at={0,0};b.at={1,0};c.at={0,1};h::hex_patch star;star.insert(a);star.insert(b);star.insert(c);++tried;
            try {
                auto tri=h::substitute_triangles(h::to_triangles(star),1);++geometric;
                h::hex_patch one;one.insert(a);auto phase=h::substitute(one,1).placement;
                (void)partial(tri,phase);accepted.try_emplace(key(star),star);
            } catch(const std::logic_error&) {}
        }
        for(unsigned steps=2;steps<=3;++steps) {
            std::map<Key,h::hex_patch> next;
            for(auto& [sig,star]:accepted) try {
                auto tri=h::substitute_triangles(h::to_triangles(star),steps);
                h::hex_patch one;one.insert(star.tiles.begin()->second);
                (void)partial(tri,h::substitute(one,steps).placement);next.emplace(sig,star);
            } catch(const std::logic_error&) {}
            accepted=std::move(next);std::cout<<"admissible_steps,"<<steps<<",types,"<<accepted.size()<<std::endl;
        }
        size_t extra=0;for(auto& [k,p]:accepted) if(!known.contains(k)) ++extra;
        std::cout<<"tried,"<<tried<<",geometrically_compatible,"<<geometric<<",accepted_types,"<<accepted.size()<<",sampled_legal,"<<known.size()<<",extra,"<<extra<<'\n';
        for(auto& [k,p]:accepted) if(!known.contains(k)) {std::cout<<"extra";for(auto [q,r,t,o]:k) std::cout<<';'<<q<<','<<r<<','<<t<<','<<o;std::cout<<'\n';}
        for(unsigned round=1;round<=5;++round) {
            std::map<Key,h::hex_patch> image;
            for(auto& [sig,star]:accepted) {
                auto tri=h::substitute_triangles(h::to_triangles(star),1);
                auto hex=partial(tri,h::substitute(star,1).placement);
                for(auto& [v,cs]:vertices(hex)) if(cs.size()==3) {
                    h::hex_patch p;
                    for(auto at:cs) {auto x=hex.at(at);x.at=at-cs.front();p.insert(x);}
                    auto sig2=key(p);if(accepted.contains(sig2)) image.try_emplace(sig2,p);
                }
            }
            accepted=std::move(image);extra=0;
            for(auto& [sig,p]:accepted) if(!known.contains(sig)) ++extra;
            std::cout<<"backwards_prune_round,"<<round<<",types,"<<accepted.size()<<",extra,"<<extra<<std::endl;
        }
        for(auto& [sig,p]:accepted) if(!known.contains(sig)) {
            std::cout<<"surviving_extra";for(auto [q,r,t,o]:sig) std::cout<<';'<<q<<','<<r<<','<<t<<','<<o;std::cout<<'\n';
        }
        auto b1=compatible_balls1(accepted);auto observed=sample_balls(depth,1);
        extra=0;for(auto& [sig,p]:b1) if(!observed.contains(sig)) ++extra;
        std::cout<<"compatible_radius1_types,"<<b1.size()<<",sampled,"<<observed.size()<<",extra,"<<extra<<std::endl;
        std::map<Key,h::hex_patch> extendible_vertices;
        for(auto& [sig,p]:b1) for(auto& [v,cs]:vertices(p)) if(cs.size()==3) {
            h::hex_patch star;for(auto at:cs) star.insert(p.at(at));extendible_vertices.try_emplace(key(star),star);
        }
        std::cout<<"extendible_vertex_types,"<<extendible_vertices.size()<<std::endl;
        std::map<Key,h::hex_patch> interior_image;
        for(auto& [sig,p]:b1) {
            auto src=h::to_triangles(p);auto moved=deform(src,g::next_lengths({}));Polygon poly;
            for(size_t j=0;j<6;++j) poly[j]=moved.at(unit[j]);
            auto tri=h::substitute_triangles(src,1);
            auto hex=partial(tri,h::substitute(p,1).placement);
            for(auto& [v,cs]:vertices(hex)) if(cs.size()==3 && inside3({3*v.q,3*v.r},poly)) {
                h::hex_patch star;for(auto at:cs) star.insert(hex.at(at));interior_image.try_emplace(key(star),star);
            }
        }
        extra=0;for(auto& [sig,p]:interior_image) if(!known.contains(sig)) ++extra;
        std::cout<<"central_face_vertex_image,"<<interior_image.size()<<",extra,"<<extra<<std::endl;
        std::set<Key> edge_types;
        for(auto& [sig,star]:interior_image) {
            for(auto a=star.tiles.begin();a!=star.tiles.end();++a) for(auto b=std::next(a);b!=star.tiles.end();++b) {
                h::hex_patch pair;pair.insert(a->second);pair.insert(b->second);edge_types.insert(key(pair));
            }
            auto tri=h::substitute_triangles(h::to_triangles(star),1);
            auto image=partial(tri,h::substitute(star,1).placement);
            for(auto& [v,cs]:vertices(image)) if(cs.size()==3) {
                h::hex_patch out;for(auto at:cs) out.insert(image.at(at));
                if(!interior_image.contains(key(out))) throw std::logic_error("certified V is not closed under substitution");
            }
        }
        std::cout<<"certified_vertex_catalog_closed,edge_types,"<<edge_types.size()<<std::endl;
        auto candidate1=compatible_balls1(interior_image);auto candidate2=compatible_balls2(candidate1);
        int universal1=999999;for(auto& [sig,p]:candidate2) universal1=std::min(universal1,geometry(p,1).radius);
        std::cout<<"local_radius1_candidates,"<<candidate1.size()<<",radius2_candidates,"<<candidate2.size()<<",min_R1,"<<universal1<<std::endl;
        std::map<Key,h::hex_patch> exact1;
        for(auto& [sig,star]:interior_image) {
            auto tri=h::substitute_triangles(h::to_triangles(star),1);
            auto hex=partial(tri,h::substitute(star,1).placement);
            for(auto& [at,x]:hex.tiles) {
                h::hex_patch p;for(auto d:ball(1)) {auto y=hex.find(at+d);if(!y) break;auto z=*y;z.at=d;p.insert(z);}
                if(p.size()==7) exact1.try_emplace(key(p),p);
            }
        }
        auto candidate2b=compatible_balls2(exact1);int universal2=999999;
        for(auto& [sig,p]:candidate2b) universal2=std::min(universal2,geometry(p,2).radius);
        std::map<Key,h::hex_patch> exact2;
        for(auto& [sig,star]:interior_image) {
            auto tri=h::substitute_triangles(h::to_triangles(star),2);
            auto hex=partial(tri,h::substitute(star,2).placement);
            for(auto& [at,x]:hex.tiles) {
                h::hex_patch p;for(auto d:ball(2)) {auto y=hex.find(at+d);if(!y) break;auto z=*y;z.at=d;p.insert(z);}
                if(p.size()==19) exact2.try_emplace(key(p),p);
            }
        }
        auto observed2=sample_balls(depth,2);extra=0;
        for(auto& [sig,p]:exact2) if(!observed2.contains(sig)) ++extra;
        std::cout<<"certified_radius1,"<<exact1.size()<<",radius2_candidates,"<<candidate2b.size()<<",min_R2,"<<universal2
                 <<",certified_radius2,"<<exact2.size()<<",radius2_extra,"<<extra<<std::endl;
        if(interior_image.size()!=19 || extra || universal1<1 || universal2<2 || exact1.size()!=14 || exact2.size()!=29)
            throw std::logic_error("certificate cardinalities changed");
        return 0;
    }
    if(mode=="catalog") {
        auto v=sample_vertices(depth);unsigned maxk=argc>3?unsigned(std::stoul(argv[3])):3;
        std::vector<std::array<std::set<Key>,6>> history;
        for(unsigned k=0;k<=maxk;++k) {
            std::array<std::set<Key>,6> catalogs;size_t minsize=999999,maxsize=0;
            for(auto& [sig,star]:v) {
                auto tri=h::substitute_triangles(h::to_triangles(star),k);
                auto phase=h::substitute(star,k).placement;
                auto hex=partial(tri,phase);minsize=std::min(minsize,hex.size());maxsize=std::max(maxsize,hex.size());
                for(int n=1;n<=6;++n) {auto b=ball(n);for(auto& [at,x]:hex.tiles) {
                    h::hex_patch collar;
                    for(auto delta:b) {auto y=hex.find(at+delta);if(!y) break;auto z=*y;z.at=delta;collar.insert(z);}
                    if(collar.size()==b.size()) catalogs[n-1].insert(key(collar));
                }}
            }
            std::cout<<"catalog,k,"<<k<<",V,"<<v.size()<<",min_hexes,"<<minsize<<",max_hexes,"<<maxsize;
            for(auto& c:catalogs) std::cout<<','<<c.size();std::cout<<std::endl;
            history.push_back(std::move(catalogs));
        }
        for(unsigned k=0;k<history.size();++k) {
            std::cout<<"equals_final,k,"<<k;
            for(unsigned n=0;n<6;++n) {
                if(maxk>=4 && !std::includes(history.back()[n].begin(),history.back()[n].end(),history[k][n].begin(),history[k][n].end()))
                    throw std::logic_error("earlier catalog not contained in certified final catalog");
                std::cout<<','<<(history[k][n]==history.back()[n]);
            }
            std::cout<<'\n';
        }
        return 0;
    }
    if(mode=="geometry") {
        auto samples=sample_balls(depth,2);unsigned maxk=argc>3?unsigned(std::stoul(argv[3])):3;
        std::cout<<"source_radius2_types,"<<samples.size()<<std::endl;
        for(unsigned k=0;k<=maxk;++k) {
            int worst=999999,best=0;Key worstkey;h::hex_patch worstsource;GeometryResult witness;
            size_t i=0,wi=0;
            for(auto& [sig,source]:samples) {
                auto r=geometry(source,k);
                if(r.radius<worst) {worst=r.radius;worstkey=sig;worstsource=source;witness=r;wi=i;}
                best=std::max(best,r.radius);++i;
            }
            std::cout<<"coverage,k,"<<k<<",min_R,"<<worst<<",max_R,"<<best<<",case,"<<wi
                     <<",target,"<<witness.target.q<<','<<witness.target.r<<",eligible_stars,"<<witness.missing.size()<<'\n';
            std::cout<<"source";for(auto [q,r,t,o]:worstkey) std::cout<<';'<<q<<','<<r<<','<<t<<','<<o;std::cout<<'\n';
            std::cout<<"actual_source";for(auto [at,x]:worstsource.tiles) std::cout<<';'<<at.q<<','<<at.r<<','<<int(x.type)<<','<<x.orientation;std::cout<<'\n';
            std::cout<<"output_phase,"<<witness.phase.origin.q<<','<<witness.phase.origin.r<<";central_polygon";
            for(auto p:witness.central) std::cout<<';'<<p.q<<','<<p.r;std::cout<<'\n';
            for(auto [v,m]:witness.missing) std::cout<<"missing,vertex,"<<v.q<<','<<v.r<<",output,"<<m.q<<','<<m.r<<'\n';
            std::cout<<std::flush;
        }
        return 0;
    }
    if(mode=="regular-red") {
        h::hex_patch p;for(auto at:ball(2)) p.insert({at,g::tile_type::red,0});
        for(unsigned k=0;k<=4;++k) {auto r=geometry(p,k);std::cout<<"regular_red,k,"<<k<<",R,"<<r.radius<<std::endl;}
        return 0;
    }
    throw std::invalid_argument("unknown probe mode");
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
