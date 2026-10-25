#include "combinatorial_hex.hpp"
#include "hex.hpp" // reference implementation is used only in this test

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

namespace g = jw::grid;
void require(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
jw::point project(g::vertex v) {
    return {static_cast<double>(v.q)+static_cast<double>(v.r)/2,
            std::sqrt(3.0)/2*static_cast<double>(v.r)};
}
void compare(const g::affine& a,const jw::matrix& b) {
    for (g::vertex v : {g::vertex{},g::vertex{1,0},g::vertex{0,1}}) {
        const auto x=project(a(v)), y=jw::apply_matrix(b,project(v));
        require(std::hypot(x.x-y.x,x.y-y.y)<1e-9,"geometric transform mismatch");
    }
}
std::set<g::cell> region(g::tile_type type,g::lengths dims,const g::affine& t={}) {
    std::set<g::cell> result;
    for (auto c:g::canonical_region(type,dims)) require(result.insert(t(c)).second,"duplicate canonical cell");
    return result;
}
int main() {
    try {
        const std::array<std::string,5> sequences{
            "blue purple yellow yellow yellow ",
            "yellow yellow yellow blue purple purple blue blue yellow yellow red ",
            "yellow yellow yellow green green purple purple purple yellow yellow yellow ",
            "yellow purple purple purple yellow yellow yellow green ",
            "yellow purple "};
        const std::array<std::size_t,5> counts{9,21,21,16,4};
        auto dims=g::lengths{};
        auto floating=jw::regular_lengths();
        for (unsigned level=0;level<4;++level) {
            const auto exact=g::rules(dims);
            const auto reference=jw::get_hex_rules(floating);
            for (auto type:g::types) {
                const auto idx=static_cast<std::size_t>(type);
                const auto& children=exact[idx];
                const auto& other=reference.at(std::string(g::name(type))).inflate_substitution;
                std::string sequence;
                std::set<g::cell> covered;
                require(children.size()==other.size(),"rule child count");
                for (std::size_t i=0;i<children.size();++i) {
                    sequence+=std::string(g::name(children[i].type))+" ";
                    compare(children[i].transform,other[i].trans);
                    for (auto c:region(children[i].type,dims,children[i].transform))
                        require(covered.insert(c).second,"rule overlap");
                }
                require(sequence==sequences[idx],"ordered rule types");
                require(covered==region(type,g::next_lengths(dims)),"rule coverage has gaps or excess");
                if (level==0) require(covered.size()==counts[idx],"depth-one triangle count");
            }
            dims=g::next_lengths(dims);
            floating=jw::next_lengths(floating);
        }
        for (auto type:g::types) for (unsigned depth : {0u,1u,2u,4u}) {
            const auto p=g::make_supertile(type,depth);
            const auto reference=jw::make_supertile(g::name(type),depth,jw::regular_lengths());
            require(p.leaves.size()==reference.size(),"recursive leaf count");
            std::set<g::cell> actual;
            for (const auto& [c,owner]:p.occupancy) {
                actual.insert(c);
                require(p.nodes.at(owner).level==0,"cell owner must be a leaf");
            }
            g::lengths outer;
            for (unsigned k=0;k<depth;++k) outer=g::next_lengths(outer);
            require(actual==region(type,outer),"recursive coverage");
            std::size_t cell_count=0;
            for (std::size_t i=0;i<p.leaves.size();++i) {
                const auto id=p.leaves[i];
                const auto& n=p.nodes.at(id);
                require(g::name(n.type)==reference[i].type,"ordered leaf type");
                compare(n.transform,reference[i].trans);
                require(std::set<g::cell>(n.cells.begin(),n.cells.end())==region(n.type,{},n.transform),"leaf cells");
                for (auto c:n.cells) require(p.occupancy.at(c)==id,"ownership mismatch");
                cell_count+=n.cells.size();
            }
            require(cell_count==p.occupancy.size(),"multiply owned cells");
            require(p.nodes[0].parent==g::no_parent && p.nodes[0].level==depth,"root identity");
            for (g::node_id id=0;id<p.nodes.size();++id) {
                const auto& n=p.nodes[id];
                if (!n.level) { require(n.children.empty(),"leaf children"); continue; }
                require(n.cells.empty(),"non-leaf owns cells");
                g::lengths d;
                for (unsigned k=1;k<n.level;++k) d=g::next_lengths(d);
                const auto rule=g::rules(d)[static_cast<std::size_t>(n.type)];
                require(n.children.size()==rule.size(),"hierarchy child count");
                for (std::size_t slot=0;slot<n.children.size();++slot) {
                    const auto child_id=n.children[slot];
                    const auto& c=p.nodes.at(child_id);
                    require(child_id>id && c.parent==id && c.child_slot==slot && c.level+1==n.level,"parent/child linkage");
                    require(c.type==rule[slot].type && c.transform==n.transform*rule[slot].transform,"ancestral frame");
                }
            }
            std::size_t exterior=0,incidences=0;
            for (const auto& [edge,items]:p.edges()) {
                require(items.size()==1 || items.size()==2,"edge multiplicity");
                incidences+=items.size();
                if (items.size()==1) ++exterior;
                for (auto item:items) {
                    require(p.occupancy.at(item.triangle)==item.owner,"edge owner");
                    const auto es=item.triangle.edges();
                    require(std::find(es.begin(),es.end(),edge)!=es.end(),"edge incidence");
                }
            }
            require(exterior>0 && incidences==3*p.occupancy.size(),"edge index completeness");
        }
        bool limited=false;
        try { (void)g::make_supertile(g::tile_type::red,std::numeric_limits<unsigned>::max()); }
        catch (const std::length_error&) { limited=true; }
        require(limited,"huge depth rejected before allocation/recursion");
        bool overflow=false;
        try { (void)g::next_lengths({std::numeric_limits<g::integer>::max(),1}); }
        catch (const std::overflow_error&) { overflow=true; }
        require(overflow,"dimension overflow");
        std::cout << "Rules, exact coverage, hierarchy and geometric comparisons passed\n";
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
