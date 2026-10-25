#include "triangle_grid.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace jw::grid {
    namespace {
        constexpr integer lo = std::numeric_limits<integer>::min();
        constexpr integer hi = std::numeric_limits<integer>::max();
        [[noreturn]] void overflow() { throw std::overflow_error("triangle-grid integer overflow"); }
        integer dot(integer a, integer b, integer c, integer d) {
            return add(multiply(a,b), multiply(c,d));
        }
    }
    integer add(integer a, integer b) {
        if ((b>0 && a>hi-b) || (b<0 && a<lo-b)) overflow();
        return a+b;
    }
    integer subtract(integer a, integer b) {
        if ((b>0 && a<lo+b) || (b<0 && a>hi+b)) overflow();
        return a-b;
    }
    integer multiply(integer a, integer b) {
        if (a>0 ? (b>0 ? a>hi/b : b<lo/a)
                : (a<0 && (b>0 ? a<lo/b : b<0 && a<hi/b))) overflow();
        return a*b;
    }
    vertex operator+(vertex a,vertex b) { return {add(a.q,b.q),add(a.r,b.r)}; }
    vertex operator-(vertex a,vertex b) { return {subtract(a.q,b.q),subtract(a.r,b.r)}; }
    edge::edge(vertex u,vertex v):a(std::min(u,v)),b(std::max(u,v)) {}
    std::array<vertex,3> cell::vertices() const {
        const auto right = anchor+vertex{1,0}, top = anchor+vertex{0,1};
        if (down) return {right,anchor+vertex{1,1},top};
        return {anchor,right,top};
    }
    std::array<edge,3> cell::edges() const {
        const auto v=vertices();
        return {edge{v[0],v[1]},edge{v[1],v[2]},edge{v[2],v[0]}};
    }
    std::array<cell,3> cell::neighbors() const {
        const integer step = down ? 1 : -1;
        return {cell{anchor,!down}, cell{anchor+vertex{step,0},!down},
                cell{anchor+vertex{0,step},!down}};
    }
    vertex affine::operator()(vertex p) const {
        return vertex{dot(a,p.q,b,p.r),dot(c,p.q,d,p.r)} + translation;
    }
    cell affine::operator()(cell t) const {
        auto vertices=t.vertices();
        for (auto& p:vertices) p=(*this)(p);
        vertex anchor=vertices[0];
        for (auto p:vertices) { anchor.q=std::min(anchor.q,p.q); anchor.r=std::min(anchor.r,p.r); }
        std::sort(vertices.begin(),vertices.end());
        for (bool down : {false,true}) {
            const cell candidate{anchor,down};
            auto v=candidate.vertices();
            std::sort(v.begin(),v.end());
            if (v==vertices) return candidate;
        }
        throw std::invalid_argument("affine map does not preserve elementary triangles");
    }
    integer affine::determinant() const { return subtract(multiply(a,d),multiply(b,c)); }
    affine affine::inverse() const {
        const auto det=determinant();
        if (det!=1 && det!=-1) throw std::invalid_argument("integer affine inverse requires determinant +/-1");
        affine result{multiply(d,det),multiply(b,-det),multiply(c,-det),multiply(a,det),{}};
        result.translation = vertex{}-result(translation);
        return result;
    }
    affine operator*(const affine& x,const affine& y) {
        return {dot(x.a,y.a,x.b,y.c), dot(x.a,y.b,x.b,y.d),
                dot(x.c,y.a,x.d,y.c), dot(x.c,y.b,x.d,y.d),x(y.translation)};
    }
    affine translate(vertex offset) { return {1,0,0,1,offset}; }
    affine rotate(int sixth_turns) {
        // R(q,r)=(-r,q+r), so R^6=I. Reduce before negating (INT_MIN is valid).
        const int k=(sixth_turns%6+6)%6;
        affine result;
        for (int i=0;i<k;++i) result=affine{0,-1,1,1,{}}*result;
        return result;
    }
    affine reflect() { return {-1,-1,0,1,{}}; }
}
