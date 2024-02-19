/*
 * IMP. Image and Mesh Processing library.
 * Copyright (C) 2016  <Konstantinos Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "IMP/IMP"


int main() {

    try {

        std::cout << "Test on 3D AHF mesh...\n";

        // Append elements
        IMP::Mesh<3,4> bs_mesh;
        bs_mesh.AppendCell({0,1,2,5});
        bs_mesh.AppendCell({0,2,5,7});
        bs_mesh.AppendCell({0,4,5,7});
        bs_mesh.AppendCell({2,5,6,7});
        bs_mesh.AppendCell({0,2,3,7});

        // now add the vertex point coordinates (5 vertices)
        bs_mesh.AppendNode({0., 0., 0.});
        bs_mesh.AppendNode({1., 0., 0.});
        bs_mesh.AppendNode({1., 1., 0.});
        bs_mesh.AppendNode({0., 1., 0.});
        bs_mesh.AppendNode({0., 0., 1.});
        bs_mesh.AppendNode({1., 0., 1.});
        bs_mesh.AppendNode({1., 1., 1.});
        bs_mesh.AppendNode({0., 1., 1.});

        bs_mesh.FindSiblingHalfFacets();
        bs_mesh.MapNodesToHalfFacets();

        // Show vertices
        std::cout << "Mesh vertices coords:\n";
        for (std::size_t i = 0; i != bs_mesh.Nodes().size(); ++i) {
            std::cout << i << " | (" << bs_mesh.Nodes()[i][0] << ", " << bs_mesh.Nodes()[i][1] << ", " << bs_mesh.Nodes()[i][2] << ")\n";
        }


        // Show elements
        std::cout << "Connectivity of all cells:\n";
        std::cout << "Cell #  |  Connectivity  |  Sibling half-facets\n";
        for (std::size_t i = 0; i != bs_mesh.Cells().size(); ++i) {
            std::cout << "   " << i << "    |      " << bs_mesh.Cells()[i] << "    |  ";
                        for (std::size_t j = 0; j != bs_mesh.Cells()[i].SibHfacets().size(); ++j) {
                            std::cout << bs_mesh.Cells()[i].SibHfacets()[j]  << " ";
                        }
                        std::cout << std::endl;
        }

        // Show vertex to half facets map.
        std::cout << "Vertex to half facets map:\n";
        for (const auto &hf : bs_mesh.NodesToHalfFacets().HalfFacets()) {
            std::cout << hf << std::endl;
        }

        // Show the connectivity of a half facet.
        std::cout << "Connectivity of the first half facet of the first cell:\n";
        std::cout << bs_mesh.Cells()[0].FacetConnectivity(0) << std::endl;

        IMP::Mesh<3,4> inp_mesh2;
        inp_mesh2.LoadFrom("/home/mood/Dropbox/a_CME_CMAME_2018/code/mesh_data/inp/cube_lvl2.inp");
        std::cout << "\n=================\n";
        IMP::Mesh<3,4> inp_mesh3;
        inp_mesh3.LoadFrom("/home/mood/Dropbox/a_CME_CMAME_2018/code/mesh_data/inp/cube_lvl3.inp");
        std::cout << "\n=================\n";
        IMP::Mesh<3,4> inp_mesh4;
        inp_mesh4.LoadFrom("/home/mood/Dropbox/a_CME_CMAME_2018/code/mesh_data/inp/cube_lvl4.inp");
        std::cout << "\n=================\n";
        IMP::Mesh<3,4> inp_mesh5;
        inp_mesh5.LoadFrom("/home/mood/Dropbox/a_CME_CMAME_2018/code/mesh_data/inp/cube_lvl5.inp");
        std::cout << "\n=================\n";
        IMP::Mesh<3,4> inp_mesh6;
        inp_mesh5.LoadFrom("/home/mood/Dropbox/a_CME_CMAME_2018/code/mesh_data/inp/cube_lvl6.inp");
        std::cout << "\n=================\n";

    }
    catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument error: " << e.what() << std::endl;
    }
    catch (const std::runtime_error &e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
    }
    catch (const std::out_of_range &e) {
        std::cerr << "Out of Range error: " << e.what() << std::endl;
    }
    catch (const std::bad_alloc &e) {
        std::cerr << "Bad allocation error:" << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception..." << std::endl;
    }

    return 0;
}
