#include "combinatorial_hex.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

namespace jw::grid {
    namespace {
        constexpr std::array<std::string_view,5> names{"yellow","blue","green","purple","red"};
        constexpr std::array<std::string_view,5> colors{"#ffff00","#3fbcec","#3da347","#6327ac","#ff0000"};
        std::size_t index(tile_type type) {
            const auto i=static_cast<std::size_t>(type);
            if (i>=types.size()) throw std::invalid_argument("unknown combinatorial tile type");
            return i;
        }
        void validate(lengths x) {
            if (!(x.L>x.S && x.S>0)) throw std::invalid_argument("integer lengths must satisfy L > S > 0");
        }
        child place(tile_type type,vertex from,vertex to,int turns,bool mirror=false) {
            // As in the source rules, 'from' is ALREADY in the reflected frame.
            const auto rotation=rotate(turns);
            return {type,translate(to-rotation(from))*rotation*(mirror ? reflect() : affine{})};
        }
    }
    std::string_view name(tile_type type) { return names.at(index(type)); }
    std::string_view color(tile_type type) { return colors.at(index(type)); }
    tile_type parse_type(std::string_view label) {
        for (auto type:types) if (name(type)==label) return type;
        throw std::invalid_argument("unknown seed type: "+std::string(label));
    }
    lengths next_lengths(lengths x) {
        validate(x);
        return {add(x.L,add(x.L,x.middle())),x.L};
    }
    rule_set rules(lengths x) {
        const auto parent=next_lengths(x);
        const integer L=x.L,M=x.middle(),A=parent.middle(),D=parent.L;
        // Symbolic basis change from the geometric rule anchors:
        // (x,y)=(q+r/2,h*r). For example (M+L/2,h*L) -> (M,L).
        // No floating-point geometry, snapping, or geometric rule calls occur here.
        const vertex O{},yellow_C{M,L},blue_C{0,M},blue_Q{L,M};
        const vertex blue_E{add(L,M),0},blue_U{subtract(D,A),A},blue_V{D,0};
        const vertex green_C1{M,0},green_U{0,A},green_P{A,0},green_C2{subtract(D,L),L};
        const vertex purple_P{L,0},purple_C{subtract(subtract(D,L),M),L},purple_V{D,0};
        const vertex red_T{0,L};
        using enum tile_type;
        rule_set result;
        result[index(yellow)]={
            place(blue,{L,0},O,-2), place(purple,O,{M,0},0),
            place(yellow,O,yellow_C,2), place(yellow,O,yellow_C,-2),
            place(yellow,O,yellow_C,1,true)};
        result[index(blue)]={
            place(yellow,O,blue_C,0), place(yellow,O,blue_C,-2), place(yellow,O,blue_C,3,true),
            place(blue,O,{M,0},0), place(purple,O,blue_Q,2), place(purple,O,blue_Q,5,true),
            place(blue,{-L,0},blue_E,1,true), place(blue,{L,0},blue_E,-2),
            place(yellow,O,blue_U,-2), place(yellow,O,blue_U,1,true), place(red,O,blue_V,2)};
        result[index(green)]={
            place(yellow,O,green_C1,2), place(yellow,O,green_C1,-1,true), place(yellow,O,green_C1,0),
            place(green,O,green_U,-2), place(green,O,green_U,3,true),
            place(purple,O,green_P,0), place(purple,O,green_P,2), place(purple,O,green_P,-1,true),
            place(yellow,O,green_C2,2), place(yellow,O,green_C2,-2), place(yellow,O,green_C2,1,true)};
        result[index(purple)]={
            place(yellow,O,O,0), place(purple,O,purple_P,0), place(purple,O,purple_P,2),
            place(purple,O,purple_P,-1,true), place(yellow,O,purple_C,-2),
            place(yellow,O,purple_C,1,true), place(yellow,O,purple_C,3,true), place(green,O,purple_V,2)};
        result[index(red)]={place(yellow,O,red_T,-2),place(purple,O,O,0)};
        return result;
    }
    std::vector<cell> canonical_region(tile_type type,lengths x) {
        validate(x);
        (void)index(type);
        const bool triangle=type==tile_type::yellow || type==tile_type::red;
        const integer height=(type==tile_type::red || type==tile_type::purple) ? x.S : x.middle();
        const integer base=triangle ? height : x.L;
        const integer count=multiply(height,add(base,subtract(base,height)));
        if (count>static_cast<integer>(max_patch_cells)) throw std::length_error("canonical region exceeds the 1,000,000-cell limit");
        std::vector<cell> result;
        result.reserve(static_cast<std::size_t>(count));
        // This row enumeration is the region definition, not polygon rasterization.
        for (integer r=0;r<height;++r) {
            const integer width=subtract(base,r);
            for (integer q=0;q<width;++q) {
                result.push_back({{q,r},false});
                if (q<width-1) result.push_back({{q,r},true});
            }
        }
        return result;
    }
    edge_index patch::edges() const {
        edge_index result;
        for (const auto& [triangle,owner]:occupancy)
            for (auto e:triangle.edges()) result[e].push_back({triangle,owner});
        return result;
    }
    patch make_supertile(tile_type seed,unsigned depth) {
        const auto seed_index=index(seed);
        const auto initial_rules=rules({});
        // Count cells before allocating or recursing. Saturation avoids count
        // overflow and rejects even UINT_MAX depth after a handful of iterations.
        std::array<std::size_t,5> counts{1,3,3,3,1};
        for (unsigned level=0;level<depth;++level) {
            std::array<std::size_t,5> next{};
            for (std::size_t i=0;i<types.size();++i)
                for (const auto& c:initial_rules[i])
                    next[i]=std::min(max_patch_cells+1,next[i]+counts[index(c.type)]);
            counts=next;
            if (counts[seed_index]>max_patch_cells)
                throw std::length_error("requested depth exceeds the 1,000,000-cell combinatorial limit");
        }
        std::vector<rule_set> levels;
        lengths dims;
        for (unsigned level=0;level<depth;++level) {
            levels.push_back(rules(dims));
            dims=next_lengths(dims);
        }
        std::array<std::vector<cell>,5> leaf_cells;
        for (auto type:types) leaf_cells[index(type)]=canonical_region(type);
        patch result;
        const auto visit=[&](auto&& self,tile_type type,unsigned level,node_id parent,
                             std::size_t slot,const affine& transform)->node_id {
            const node_id id=result.nodes.size();
            result.nodes.push_back({type,level,parent,slot,transform,{},{}});
            if (level==0) {
                result.leaves.push_back(id);
                for (auto local:leaf_cells[index(type)]) {
                    const auto world=transform(local);
                    if (!result.occupancy.emplace(world,id).second)
                        throw std::logic_error("combinatorial substitution produced overlapping cells");
                    result.nodes[id].cells.push_back(world);
                }
            } else {
                const auto& children=levels[level-1][index(type)];
                for (std::size_t i=0;i<children.size();++i) {
                    const auto& child=children[i];
                    const auto child_id=self(self,child.type,level-1,id,i,transform*child.transform);
                    // Recursive appends can reallocate nodes; retain IDs, not references.
                    result.nodes[id].children.push_back(child_id);
                }
            }
            return id;
        };
        visit(visit,seed,depth,no_parent,0,affine{});
        if (result.occupancy.size()!=counts[seed_index]) throw std::logic_error("incorrect substitution cell count");
        return result;
    }
}
