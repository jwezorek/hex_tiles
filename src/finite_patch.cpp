#include "finite_patch.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <queue>
#include <set>
#include <optional>
#include <stdexcept>
#include <vector>

namespace jw::hexgrid {
    namespace {
        using grid::affine;
        using grid::cell;
        using grid::integer;
        using grid::node_id;
        using grid::occurrence;
        using grid::patch;
        using grid::tile_type;
        using grid::vertex;

        // Each canonical polygon is specified by its corner roles. The
        // matching corner in the next length realization is the same role,
        // NOT the same lattice point or a uniformly scaled point.
        std::vector<vertex> corners(tile_type type,grid::lengths d) {
            const auto i=static_cast<unsigned>(type);
            if (i>=grid::types.size()) throw std::invalid_argument("invalid triangle tile type");
            const bool small=type==tile_type::yellow || type==tile_type::red;
            const integer h=(type==tile_type::red || type==tile_type::purple)?d.S:d.middle();
            if (small) return {{0,0},{h,0},{0,h}};
            return {{0,0},{d.L,0},{grid::subtract(d.L,h),h},{0,h}};
        }
        // Index the pairs of source leaves sharing an elementary edge. Used
        // to check that new placement did not open a former interface.
        using pairs=std::set<std::pair<node_id,node_id>>;
        pairs shared_edges(const std::map<cell,node_id>& occupancy) {
            std::map<grid::edge,node_id> first;
            pairs result;
            for (const auto& [c,id]:occupancy) for (auto e:c.edges()) {
                const auto [it,inserted]=first.emplace(e,id);
                if (!inserted && it->second!=id)
                    result.emplace(std::min(id,it->second),std::max(id,it->second));
            }
            return result;
        }
        std::int64_t euler_characteristic(const std::map<cell,node_id>& occupancy) {
            std::vector<vertex> vertices;
            std::vector<grid::edge> edges;
            vertices.reserve(occupancy.size()*3);
            edges.reserve(occupancy.size()*3);
            for (const auto& [triangle,id]:occupancy) {
                (void)id;
                for (auto p:triangle.vertices()) vertices.push_back(p);
                for (auto e:triangle.edges()) edges.push_back(e);
            }
            std::sort(vertices.begin(),vertices.end());
            vertices.erase(std::unique(vertices.begin(),vertices.end()),vertices.end());
            std::sort(edges.begin(),edges.end());
            edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
            // V - E + F is invariant under a valid deformation of a
            // finite union of triangular cells. For a connected patch it
            // detects the creation or elimination of unintended holes.
            return static_cast<std::int64_t>(vertices.size())-
                   static_cast<std::int64_t>(edges.size())+
                   static_cast<std::int64_t>(occupancy.size());
        }
        void validate(const patch& input) {
            if (input.occupancy.size()>grid::max_patch_cells)
                throw std::length_error("triangle patch exceeds cell limit");
            if (input.leaves.size()!=input.nodes.size())
                throw std::invalid_argument("finite patch must consist of base-length leaves");
            std::map<cell,node_id> verified;
            for (auto id:input.leaves) {
                if (id>=input.nodes.size()) throw std::invalid_argument("invalid leaf ID");
                const auto& n=input.nodes[id];
                if (n.level!=0 ||
                    (n.transform.determinant()!=1 && n.transform.determinant()!=-1))
                    throw std::invalid_argument("invalid finite-patch leaf transform");
                if (n.cells.size()!=grid::canonical_region(n.type).size())
                    throw std::invalid_argument("incomplete base-length prototile");
                for (auto local:grid::canonical_region(n.type))
                    if (!verified.emplace(n.transform(local),id).second)
                        throw std::invalid_argument("overlapping base-length prototiles");
            }
            if (verified!=input.occupancy)
                throw std::invalid_argument("triangle patch occupancy/leaf mismatch");
        }
        struct corner_incidence {
            node_id id;
            vertex new_local; // transformed by the old dihedral linear map only
        };
        patch once(const patch& input) {
            validate(input);
            if (input.leaves.empty()) return {};
            const auto new_dims=grid::next_lengths({});
            const auto rewrite=grid::rules({});
            std::map<vertex,std::vector<corner_incidence>> junctions;
            for (auto id:input.leaves) {
                const auto& n=input.nodes[id];
                const affine linear{n.transform.a,n.transform.b,n.transform.c,n.transform.d,{}};
                const auto old=corners(n.type,{}),updated=corners(n.type,new_dims);
                for (std::size_t k=0;k<old.size();++k)
                    junctions[n.transform(old[k])].push_back({id,linear(updated[k])});
            }
            struct constraint {node_id to;vertex delta;};
            std::vector<std::vector<constraint>> graph(input.nodes.size());
            for (const auto& [point,items]:junctions) {
                (void)point;
                if (items.size()<2) continue;
                const auto a=items[0];
                for (std::size_t k=1;k<items.size();++k) {
                    const auto b=items[k];
                    if (a.id==b.id) continue;
                    // new_offset[b] = new_offset[a] + a.local - b.local
                    const auto delta=a.new_local-b.new_local;
                    graph[a.id].push_back({b.id,delta});
                    graph[b.id].push_back({a.id,vertex{}-delta});
                }
            }
            std::vector<std::optional<vertex>> offsets(input.nodes.size());
            std::queue<node_id> pending;
            std::size_t components=0;
            // The absolute position of a component is an arbitrary gauge:
            // anchor its lexicographically first source leaf at its old origin.
            for (auto root:input.leaves) {
                if (offsets[root]) continue;
                ++components;
                offsets[root]=input.nodes[root].transform.translation;
                pending.push(root);
                while (!pending.empty()) {
                    const auto at=pending.front();pending.pop();
                    for (auto [to,delta]:graph[at]) {
                        const auto desired=*offsets[at]+delta;
                        if (!offsets[to]) {offsets[to]=desired;pending.push(to);}
                        else if (*offsets[to]!=desired)
                            throw std::invalid_argument("incompatible shared-corner constraints; patch cannot be substituted");
                    }
                }
            }
            if (components>1)
                throw std::invalid_argument("finite triangle patch is not corner-connected; relative placement is undetermined");

            std::size_t expanded_count=0;
            for (auto id:input.leaves) {
                const auto& n=input.nodes[id];
                expanded_count+=grid::canonical_region(n.type,new_dims).size();
                if (expanded_count>grid::max_patch_cells)
                    throw std::length_error("substitution exceeds 1,000,000-cell limit");
            }
            patch output;
            std::map<cell,node_id> parent_occupancy;
            for (auto id:input.leaves) {
                const auto& n=input.nodes[id];
                const affine parent{n.transform.a,n.transform.b,n.transform.c,n.transform.d,*offsets[id]};
                for (auto c:grid::canonical_region(n.type,new_dims)) {
                    if (!parent_occupancy.emplace(parent(c),id).second)
                        throw std::invalid_argument("expanded parent prototiles overlap");
                }
                const auto& children=rewrite.at(static_cast<unsigned>(n.type));
                for (std::size_t slot=0;slot<children.size();++slot) {
                    const auto& c=children[slot];
                    const auto child_transform=parent*c.transform;
                    const auto child_id=output.nodes.size();
                    output.nodes.push_back({c.type,0,grid::no_parent,slot,child_transform,{},{}});
                    output.leaves.push_back(child_id);
                    for (auto local:grid::canonical_region(c.type)) {
                        const auto world=child_transform(local);
                        if (!output.occupancy.emplace(world,child_id).second)
                            throw std::logic_error("substituted child triangles overlap");
                        output.nodes.back().cells.push_back(world);
                    }
                }
            }
            if (output.occupancy.size()!=parent_occupancy.size())
                throw std::logic_error("substitution left gaps inside expanded parents");
            for (const auto& [c,id]:parent_occupancy) {
                (void)id;
                if (!output.occupancy.contains(c))
                    throw std::logic_error("substitution did not cover expanded parent");
            }
            if (shared_edges(input.occupancy)!=shared_edges(parent_occupancy))
                throw std::invalid_argument("substitution does not preserve subtile edge adjacencies");
            if (euler_characteristic(input.occupancy)!=euler_characteristic(parent_occupancy))
                throw std::invalid_argument("substitution changes patch topology or creates a hole");
            return output;
        }
    }
    grid::patch substitute_triangles(const grid::patch& input,unsigned steps) {
        validate(input);
        if (steps==0 || input.leaves.empty()) return input;
        grid::patch result=input;
        for (unsigned step=0;step<steps;++step) result=once(result);
        return result;
    }
    hex_patch substitute(const hex_patch& input,unsigned steps) {
        const auto triangles=to_triangles(input);
        if (!steps) {
            auto result=convert(triangles,input.placement);
            if (result.size()!=input.size())
                throw std::logic_error("hex patch failed zero-step reconstruction");
            return result;
        }
        // Determine the output phase from an independently substituted copy
        // of one complete input hex. Its first leaf is also the first leaf of
        // the full (map-ordered) input patch, so the shared-corner solver fixes
        // the SAME translational gauge at every step. Exactly one of the three
        // lattice phases must admit valid complete hex motifs in this sample.
        if (input.empty()) {
            hex_patch empty;
            empty.placement=input.placement;
            return empty;
        }
        hex_patch witness;
        witness.placement=input.placement;
        witness.insert(input.tiles.begin()->second);
        const auto sample=substitute_triangles(to_triangles(witness),steps);
        std::optional<lattice> phase;
        for (integer k=0;k<3;++k) {
            lattice candidate{{k,0}};
            try {
                const auto trial=convert(sample,candidate);
                if (!trial.empty()) {
                    if (phase) throw std::invalid_argument("ambiguous substituted hex lattice phase");
                    phase=candidate;
                }
            } catch (const std::logic_error&) {
                // An off-phase ring has a complete but invalid motif.
            }
        }
        if (!phase) throw std::logic_error("substitution has no valid hex lattice phase");
        auto expanded=substitute_triangles(triangles,steps);
        return convert(expanded,*phase);
    }
}
