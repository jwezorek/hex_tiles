#include "hex_grid.hpp"

#include <algorithm>
#include <array>
#include <set>
#include <stdexcept>
#include <string>

namespace jw::hexgrid {
    using grid::integer;
    using grid::vertex;
    using grid::cell;

    coordinate operator+(coordinate a,coordinate b) {
        return {grid::add(a.q,b.q),grid::add(a.r,b.r)};
    }
    coordinate operator-(coordinate a,coordinate b) {
        return {grid::subtract(a.q,b.q),grid::subtract(a.r,b.r)};
    }
    std::array<coordinate,6> neighbors(coordinate at) {
        std::array<coordinate,6> result;
        for (std::size_t i=0;i<6;++i) result[i]=at+directions[i];
        return result;
    }
    vertex lattice::center(coordinate at) const {
        return origin+vertex{grid::subtract(at.q,at.r),
                             grid::add(at.q,grid::multiply(2,at.r))};
    }
    coordinate lattice::locate(vertex v) const {
        const auto delta=v-origin;
        const auto triple_r=grid::subtract(delta.r,delta.q);
        if (triple_r%3) throw std::invalid_argument("triangle vertex is not a hex-grid center");
        const integer r=triple_r/3;
        return {grid::add(delta.q,r),r};
    }
    bool hex_patch::insert(hex_tile tile) {
        return tiles.emplace(tile.at,std::move(tile)).second;
    }
    const hex_tile* hex_patch::find(coordinate at) const {
        const auto i=tiles.find(at);
        return i==tiles.end()?nullptr:&i->second;
    }
    const hex_tile& hex_patch::at(coordinate where) const { return tiles.at(where); }
    std::optional<std::pair<coordinate,coordinate>> hex_patch::bounds() const {
        if (tiles.empty()) return std::nullopt;
        coordinate low{tiles.begin()->first},high{low};
        for (const auto& [where,tile]:tiles) {
            (void)tile;
            low.q=std::min(low.q,where.q); low.r=std::min(low.r,where.r);
            high.q=std::max(high.q,where.q); high.r=std::max(high.r,where.r);
        }
        return std::pair{low,high};
    }

    namespace {
        // A triangle's three vertices occupy all three residues of (q-r) mod 3.
        // Therefore exactly one vertex belongs to any selected hex-center lattice:
        // the six incident triangles partition into disjoint candidate hexes.
        int residue(vertex p) {
            // Avoid subtraction overflow even at extreme valid triangle coordinates.
            const int q=static_cast<int>((p.q%3+3)%3);
            const int r=static_cast<int>((p.r%3+3)%3);
            return (q-r+3)%3;
        }
        std::array<cell,6> around(vertex p) {
            const vertex below=p+vertex{0,-1}, left=p+vertex{-1,0};
            return {{ {p+vertex{-1,-1},true}, {below,false}, {below,true},
                      {p,false}, {left,true}, {left,false} }};
        }
        struct decoration {
            int direction=0;
            bool mirrored=false;
            auto operator<=>(const decoration&) const = default;
        };
        using signature=std::array<decoration,6>;
        // Canonical six-sector leaf-arrow motifs, from the exact unit-subtile rules.
        // Triangle hexes are six 1-cell subtiles. Trapezoid hexes are two 3-cell
        // subtiles. In both cases the source's +q arrow is the decoration direction.
        constexpr signature triangle_motif{{{4,true},{4,false},{0,true},
                                             {0,false},{2,true},{2,false}}};
        constexpr signature trapezoid_motif{{{0,true},{0,true},{0,true},
                                              {0,false},{0,false},{0,false}}};
        constexpr std::array<vertex,6> arrows{{{1,0},{0,1},{-1,1},
                                                {-1,0},{0,-1},{1,-1}}};
        signature rotate_signature(const signature& from,int turns) {
            signature result;
            for (int i=0;i<6;++i) {
                const auto& d=from[i];
                result[(i+turns)%6]={ (d.direction+turns)%6, d.mirrored };
            }
            return result;
        }
        std::optional<unsigned> motif_orientation(grid::tile_type type,const signature& observed) {
            const auto& template_motif=(type==grid::tile_type::yellow || type==grid::tile_type::red)
                                      ?triangle_motif:trapezoid_motif;
            // The first matching rotation is canonical for motifs with symmetry.
            for (unsigned turn=0;turn<6;++turn)
                if (rotate_signature(template_motif,static_cast<int>(turn))==observed) return turn;
            return std::nullopt;
        }
        int direction(vertex v) {
            for (int i=0;i<6;++i) if (arrows[i]==v) return i;
            throw std::logic_error("hex reconstruction: leaf arrow is not a unit lattice direction");
        }
        struct candidate {
            std::array<cell,6> cells;
            std::array<grid::node_id,6> ids{};
            int occupied=0;
            std::optional<grid::tile_type> type;
            std::optional<unsigned> orientation;
        };
        candidate examine(const grid::patch& p,vertex center) {
            candidate c;
            c.cells=around(center);
            signature observed;
            for (int k=0;k<6;++k) {
                const auto it=p.occupancy.find(c.cells[k]);
                if (it==p.occupancy.end()) continue;
                ++c.occupied;
                const auto id=it->second;
                if (id>=p.nodes.size()) throw std::logic_error("hex reconstruction: invalid leaf ID");
                const auto& owner=p.nodes[id];
                if (owner.level!=0) throw std::logic_error("hex reconstruction: non-leaf triangle owner");
                if (!c.type) c.type=owner.type;
                if (c.type!=owner.type) c.type.reset(); // fails the motif check below
                c.ids[k]=id;
                const auto determinant=owner.transform.determinant();
                if (determinant!=1 && determinant!=-1)
                    throw std::logic_error("hex reconstruction: invalid leaf affine determinant");
                observed[k]={direction({owner.transform.a,owner.transform.c}),determinant<0};
            }
            if (c.occupied==6) {
                // The subtype must be uniform; recheck explicitly because an
                // arbitrary alternating sequence can restore an optional value.
                const auto type=p.nodes.at(c.ids[0]).type;
                const bool uniform=std::all_of(c.ids.begin(),c.ids.end(),[&](auto id){
                    return p.nodes.at(id).type==type;
                });
                if (uniform) {
                    c.type=type;
                    // Check that a unit triangle appears once per source leaf,
                    // or a three-cell trapezoid appears exactly twice.
                    std::map<grid::node_id,int> multiplicity;
                    for (auto id:c.ids) ++multiplicity[id];
                    const bool tiny=type==grid::tile_type::yellow || type==grid::tile_type::red;
                    const bool correct_owners=tiny ? multiplicity.size()==6 :
                        (multiplicity.size()==2 &&
                         std::all_of(multiplicity.begin(),multiplicity.end(),
                             [](const auto& entry){return entry.second==3;}));
                    if (correct_owners) c.orientation=motif_orientation(type,observed);
                }
            }
            return c;
        }
        // The root canonical region has an exact integer predicate. This also
        // works after any dihedral transformation of an entire patch.
        struct support {
            grid::tile_type type;
            grid::lengths dims;
            grid::affine inverse;
            explicit support(const grid::patch& p)
                :type(p.nodes.at(0).type),inverse(p.nodes.at(0).transform.inverse()) {
                for (unsigned n=0;n<p.nodes.at(0).level;++n) dims=grid::next_lengths(dims);
            }
            bool contains(cell world) const {
                const auto local=inverse(world);
                const bool tiny=type==grid::tile_type::yellow || type==grid::tile_type::red;
                const integer height=(type==grid::tile_type::red || type==grid::tile_type::purple)
                    ?dims.S:dims.middle();
                const integer base=tiny?height:dims.L;
                const auto q=local.anchor.q,r=local.anchor.r;
                if (q<0 || r<0 || r>=height || q>=base) return false;
                // q+r < base for up, q+r < base-1 for down.
                return q<base-r-(local.down?1:0);
            }
            bool strictly_inside(cell c) const {
                if (!contains(c)) return false;
                for (auto n:c.neighbors()) if (!contains(n)) return false;
                return true;
            }
        };
    }

    hex_patch convert(const grid::patch& source) {
        if (source.nodes.empty()) throw std::invalid_argument("cannot convert a patch without a root");
        const support expected(source);
        // Verify that the input triangle occupancy is a subset of the root
        // region, and each claimed leaf actually contains its own triangle.
        for (const auto& [c,id]:source.occupancy) {
            if (!expected.contains(c)) throw std::logic_error("triangle outside the root supertile");
            if (id>=source.nodes.size() || source.nodes[id].level!=0 ||
                std::find(source.nodes[id].cells.begin(),source.nodes[id].cells.end(),c)==source.nodes[id].cells.end())
                throw std::logic_error("triangle occupancy disagrees with leaf cells");
        }
        std::array<std::set<vertex>,3> centers;
        for (const auto& [c,id]:source.occupancy) {
            (void)id;
            for (auto v:c.vertices()) centers[residue(v)].insert(v);
        }
        // Select the global hex-grid phase by complete, decorated motifs.
        // A patch can contain a few apparently monochrome off-phase rings;
        // only the true phase forms a coherent tiling of decorated hexes.
        std::array<std::size_t,3> valid{},invalid{};
        for (int phase=0;phase<3;++phase) for (auto center:centers[phase]) {
            const auto c=examine(source,center);
            if (c.occupied==6) {
                if (c.orientation) ++valid[phase];
                else ++invalid[phase];
            }
        }
        int chosen=0;
        for (int phase=1;phase<3;++phase)
            if (valid[phase]>valid[chosen] ||
                (valid[phase]==valid[chosen] && invalid[phase]<invalid[chosen]))
                chosen=phase;
        if (valid[chosen]>0) for (int phase=0;phase<3;++phase)
            if (phase!=chosen && valid[phase]==valid[chosen] && invalid[phase]==invalid[chosen])
                throw std::logic_error("ambiguous hex-grid phase in patch");
        hex_patch result;
        result.placement.origin={chosen,0};
        std::set<cell> used;
        std::size_t boundary_cells=0;
        for (auto center:centers[chosen]) {
            auto c=examine(source,center);
            if (c.occupied!=6) {
                // A partial ring crossing the outside is expected. A missing
                // elementary cell surrounded by the ROOT's intended region is not.
                for (const auto& missing:c.cells)
                    if (!source.occupancy.contains(missing) && expected.strictly_inside(missing))
                        throw std::logic_error("incomplete interior hexagon at triangle vertex ("+
                            std::to_string(center.q)+","+std::to_string(center.r)+")");
                ++result.rejected_boundary;
                boundary_cells+=static_cast<std::size_t>(c.occupied);
                continue;
            }
            if (!c.orientation) throw std::logic_error("inconsistent complete hexagon at triangle vertex ("+
                std::to_string(center.q)+","+std::to_string(center.r)+")");
            for (const auto& triangle:c.cells)
                if (!used.insert(triangle).second)
                    throw std::logic_error("triangle assigned to multiple hexagons");
            hex_tile tile{result.placement.locate(center),*c.type,*c.orientation,c.cells,c.ids};
            if (!result.insert(std::move(tile)))
                throw std::logic_error("duplicate hexagon coordinate");
        }
        if (used.size()+boundary_cells!=source.occupancy.size())
            throw std::logic_error("triangle cells are not partitioned by the hex lattice");
        return result;
    }
}
