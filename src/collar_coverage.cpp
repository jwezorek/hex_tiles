// Exact coverage reduction from research/collar-depth.md equation (4).
// Every complete three-face star in each certified context is considered.
#include "collars_detail.hpp"
#include <algorithm>
#include <map>
#include <limits>
#include <mutex>
#include <queue>
#include <set>
#include <stdexcept>
namespace jw::hexgrid::detail {
namespace {
namespace g=jw::grid;namespace h=jw::hexgrid;
using P=g::vertex;using C=h::coordinate;
static int mod(int a,int b) {return (a%b+b)%b;}
static std::vector<P> corners(g::tile_type t,g::lengths d) {
    bool tiny=t==g::tile_type::yellow || t==g::tile_type::red;
    auto z=(t==g::tile_type::purple || t==g::tile_type::red)?d.S:d.middle();
    if(tiny) return {{0,0},{z,0},{0,z}};
    return {{0,0},{d.L,0},{g::subtract(d.L,z),z},{0,z}};
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

static std::map<P,std::vector<C>> vertices(const h::hex_patch& p) {
    std::map<P,std::vector<C>> out;
    for(auto& [at,x]:p.tiles) for(auto u:unit) out[p.placement.center(at)+u].push_back(at);
    return out;
}

using Polygon=std::array<P,6>;
static bool inside3(P triple,const Polygon& polygon) {
    for(size_t i=0;i<6;++i) {
        P a=polygon[i],d=polygon[(i+1)%6]-a;
        auto x=triple-P{g::multiply(3,a.q),g::multiply(3,a.r)};
        if(g::subtract(g::multiply(d.q,x.r),g::multiply(d.r,x.q))<0) return false;
    }
    return true;
}
static std::set<C> footprint(const Polygon& polygon,h::lattice lattice) {
    long long qmin=polygon[0].q,qmax=qmin,rmin=polygon[0].r,rmax=rmin;
    for(auto p:polygon) {qmin=std::min(qmin,p.q);qmax=std::max(qmax,p.q);rmin=std::min(rmin,p.r);rmax=std::max(rmax,p.r);}
    std::set<C> out;
    const auto lowq=g::subtract(qmin,1),highq=g::add(qmax,1);
    const auto lowr=g::subtract(rmin,1),highr=g::add(rmax,1);
    const auto estimate=g::multiply(g::add(g::subtract(highq,lowq),1),g::add(g::subtract(highr,lowr),1));
    if(estimate>4'000'000) throw std::length_error("coverage footprint exceeds scan budget");
    for(auto q=lowq;q<=highq;q=g::add(q,1)) for(auto r=lowr;r<=highr;r=g::add(r,1)) {
        P p{q,r};C at;
        if(mod(int(q%3-r%3-lattice.origin.q%3+lattice.origin.r%3),3)) continue;
        at=lattice.locate(p);
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

static g::integer context_radius(const h::hex_patch& source,unsigned k) {
    auto src=h::to_triangles(source);g::lengths dims;for(unsigned i=0;i<k;++i) dims=g::next_lengths(dims);
    auto moved=deform(src,dims);auto phase=k?direct_phase(src,moved,k):source.placement;
    std::map<C,Polygon> polygons;std::map<C,std::set<C>> images;
    for(auto& [at,tile]:source.tiles) {
        Polygon poly;for(size_t i=0;i<6;++i) poly[i]=moved.at(source.placement.center(at)+unit[i]);
        polygons[at]=poly;images[at]=footprint(poly,phase);
    }
    std::map<C,g::integer> best;
    for(auto& [v,cs]:vertices(source)) if(cs.size()==3) {
        std::set<C> image;for(auto a:cs) image.insert(images[a].begin(),images[a].end());
        std::map<C,g::integer> distance;std::queue<C> todo;
        for(auto x:image) {
            auto ns=h::neighbors(x);
            if(std::any_of(ns.begin(),ns.end(),[&](C y){return !image.contains(y);})) {distance[x]=0;todo.push(x);}
        }
        while(!todo.empty()) {auto x=todo.front();todo.pop();for(auto y:h::neighbors(x)) {
            if(image.contains(y) && !distance.contains(y)) {distance[y]=g::add(distance[x],1);todo.push(y);}
        }}
        for(auto [x,d]:distance) {
            auto found=best.find(x);if(found==best.end()) best.emplace(x,d);else found->second=std::max(found->second,d);
        }
    }
    g::integer result=std::numeric_limits<g::integer>::max();
    for(auto x:images.at({})) {
        P p=phase.center(x);
        if(!inside3({g::multiply(3,p.q),g::multiply(3,p.r)},polygons.at({}))) continue;
        result=std::min(result,best.contains(x)?best[x]:-1);
    }
    if(result==std::numeric_limits<g::integer>::max()) throw std::logic_error("coverage polygon has no output centers");
    return result;
}
}
grid::integer coverage_radius(unsigned depth) {
    if(depth>5) throw std::length_error("coverage depth exceeds supported limit 5");
    static std::mutex lock;
    static std::array<std::optional<grid::integer>,6> cache;
    const std::lock_guard guard(lock);
    if(cache[depth]) return *cache[depth];
    grid::integer result=std::numeric_limits<grid::integer>::max();
    for(const auto& source:certified_radius2_contexts()) result=std::min(result,context_radius(source,depth));
    cache[depth]=result;return result;
}
}
