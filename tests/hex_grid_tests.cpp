#include "hex_grid.hpp"
#include "hex_grid_svg.hpp"

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

namespace g = jw::grid;
namespace h = jw::hexgrid;

namespace {
    void require(bool ok,const char* message) {
        if (!ok) throw std::runtime_error(message);
    }
    template<typename F> void expect_failure(F action,const char* message) {
        try { action(); } catch (const std::logic_error&) { return; }
        throw std::runtime_error(message);
    }
    g::patch transformed(g::patch p,g::affine t) {
        for (auto& node:p.nodes) {
            node.transform=t*node.transform;
            for (auto& c:node.cells) c=t(c);
        }
        std::map<g::cell,g::node_id> moved;
        for (const auto& [c,id]:p.occupancy)
            require(moved.emplace(t(c),id).second,"transformation collision");
        p.occupancy=std::move(moved);
        return p;
    }
    void delete_cell(g::patch& p,g::cell c) {
        const auto id=p.occupancy.at(c);
        require(p.occupancy.erase(c)==1,"remove test triangle");
        auto& cells=p.nodes.at(id).cells;
        const auto it=std::find(cells.begin(),cells.end(),c);
        require(it!=cells.end(),"triangle belongs to its owner");
        cells.erase(it);
    }
    // A source triangle is strictly inside a generated patch if each of its
    // neighbors also exists. Choose a cell in a hex with six hex neighbors.
    const h::hex_tile& interior(const h::hex_patch& p) {
        for (const auto& [where,tile]:p.tiles) {
            auto n=h::neighbors(where);
            if (std::all_of(n.begin(),n.end(),[&](auto at){return p.find(at)!=nullptr;}))
                return tile;
        }
        throw std::runtime_error("test patch has no interior hex");
    }
    void validate_patch(const g::patch& source,const h::hex_patch& hp) {
        std::set<g::cell> assigned;
        require(hp.size()*6<=source.occupancy.size(),"too many reconstructed cells");
        for (const auto& [at,tile]:hp.tiles) {
            require(tile.at==at && hp.find(at)==&tile,"indexed tile identity");
            require(hp.placement.locate(hp.placement.center(at))==at,"hex lattice round trip");
            require(tile.orientation<6,"discrete orientation");
            const auto center=hp.placement.center(at);
            for (std::size_t j=0;j<6;++j) {
                const auto c=tile.triangles[j];
                const auto points=c.vertices();
                require(std::find(points.begin(),points.end(),center)!=points.end(),
                        "each sector touches the corresponding hex center");
                require(source.occupancy.at(c)==tile.source_leaves[j],"triangle source identity");
                require(source.nodes.at(tile.source_leaves[j]).type==tile.type,"triangle color identity");
                require(assigned.insert(c).second,"triangles assigned twice");
            }
        }
        require(assigned.size()==hp.size()*6,"each complete tile has six distinct triangles");
        if (hp.empty()) require(!hp.bounds(),"empty bounds");
        else {
            const auto limits=*hp.bounds();
            for (const auto& [at,tile]:hp.tiles) {
                (void)tile;
                require(at.q>=limits.first.q && at.q<=limits.second.q &&
                        at.r>=limits.first.r && at.r<=limits.second.r,"hex axial bounds");
            }
        }
    }
}

int main() {
    try {
        using h::coordinate;
        const coordinate base{8,-5};
        require(base+coordinate{-4,7}==coordinate{4,2},"hex coordinate addition");
        require(base-coordinate{-4,7}==coordinate{12,-12},"hex coordinate subtraction");
        const auto around=h::neighbors(base);
        for (std::size_t i=0;i<6;++i)
            require(around[i]==base+h::directions[i],"hex neighbor offsets");
        h::lattice lattice{{2,-3}};
        for (g::integer q=-10;q<=10;++q) for (g::integer r=-10;r<=10;++r)
            require(lattice.locate(lattice.center({q,r}))==coordinate{q,r},"exact lattice conversion");
        bool off_grid=false;
        try { (void)lattice.locate(lattice.origin+g::vertex{1,0}); }
        catch (const std::invalid_argument&) { off_grid=true; }
        require(off_grid,"off-grid vertex rejected");

        constexpr std::array<std::array<std::size_t,5>,5> expected{{
            {{0,0,7,63,477}}, {{0,2,18,150,1089}}, {{0,1,18,150,1089}},
            {{0,1,12,106,781}}, {{0,0,2,22,176}}
        }};
        for (auto type:g::types) for (unsigned depth=0;depth<=4;++depth) {
            const auto p=g::make_supertile(type,depth);
            const auto hp=h::convert(p);
            require(hp.size()==expected[static_cast<unsigned>(type)][depth],
                    "substitution-depth hexagon count");
            validate_patch(p,hp);
        }

        const auto source=g::make_supertile(g::tile_type::green,4);
        const auto original=h::convert(source);
        const auto& mid=interior(original);
        for (unsigned turns=0;turns<6;++turns) for (bool mirror:{false,true}) {
            const auto transform=g::translate({-8,13})*g::rotate(static_cast<int>(turns))*
                                  (mirror?g::reflect():g::affine{});
            const auto p=transformed(source,transform);
            const auto changed=h::convert(p);
            require(changed.size()==original.size(),"rigid symmetry preserves tile count");
            require(changed.rejected_boundary==original.rejected_boundary,
                    "rigid symmetry preserves boundary fragments");
            validate_patch(p,changed);
            for (const auto& [coord,tile]:original.tiles) {
                const auto destination=transform(original.placement.center(coord));
                const auto new_at=changed.placement.locate(destination);
                const auto& actual=changed.at(new_at);
                require(tile.type==actual.type,"rigid symmetry preserves hex type");
                const bool tiny=tile.type==g::tile_type::yellow || tile.type==g::tile_type::red;
                const int modulus=tiny?2:6;
                const int raw=mirror ? (tiny?1:3)+static_cast<int>(turns)-
                                         static_cast<int>(tile.orientation)
                                     : static_cast<int>(tile.orientation+turns);
                const unsigned predicted=static_cast<unsigned>((raw%modulus+modulus)%modulus);
                require(actual.orientation==predicted,"rigid symmetry transforms hex orientation");
                for (std::size_t j=0;j<6;++j)
                    require(p.occupancy.contains(transform(tile.triangles[j])),
                            "rigid symmetry preserves constituent triangles");
            }
        }
        // Delete a source triangle in the center of a fully surrounded hex.
        auto damaged=source;
        delete_cell(damaged,mid.triangles[0]);
        expect_failure([&]{ (void)h::convert(damaged); },"internal missing triangle not reported");

        // Deleting an edge triangle creates only another exterior fragment.
        auto clipped=source;
        bool removed=false;
        for (const auto& [c,id]:clipped.occupancy) {
            (void)id;
            const auto ns=c.neighbors();
            if (std::any_of(ns.begin(),ns.end(),[&](auto n){return !clipped.occupancy.contains(n);})) {
                delete_cell(clipped,c);
                removed=true;
                break;
            }
        }
        require(removed,"test patch perimeter");
        const auto boundary=h::convert(clipped);
        require(boundary.size()<=original.size(),"cut boundary cannot add complete hexes");
        require(boundary.rejected_boundary>0,"cut boundary fragments reported");
        validate_patch(clipped,boundary);

        // A complete interior tile containing mixed subtiles is not acceptable.
        auto mixed=source;
        const auto colored_it=std::find_if(original.tiles.begin(),original.tiles.end(),
            [](const auto& kv){return kv.second.type==g::tile_type::yellow;});
        require(colored_it!=original.tiles.end(),"find yellow tile for mixed test");
        const auto& colored=colored_it->second;
        auto& leaf=mixed.nodes.at(colored.source_leaves[0]);
        leaf.type=g::tile_type::blue;
        expect_failure([&]{ (void)h::convert(mixed); },"interior mixed hexagon not rejected");

        h::hex_patch independent;
        independent.placement=original.placement;
        require(independent.insert(mid),"first tile inserted");
        require(!independent.insert(mid),"duplicate coordinate rejected");
        require(independent.size()==1,"duplicate insertion leaves patch unchanged");

        const auto path=std::filesystem::temp_directory_path()/"hex_grid_reconstruction_test.svg";
        h::to_svg(path.string(),original);
        std::ifstream in(path);
        const std::string xml{std::istreambuf_iterator<char>(in),{}};
        require(xml.find("id=\"hexagons\"")!=std::string::npos,"hex SVG group");
        require(xml.find("data-orientation=\"")!=std::string::npos,"hex SVG metadata");
        require(xml.find("marker-end=\"url(#hex-arrowhead)\"")!=std::string::npos,"hex SVG arrows");
        std::filesystem::remove(path);
        std::cout<<"Hexagon coordinates, reconstruction, boundaries, ornaments, rotations, "
                    "reflections and SVG passed\n";
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
