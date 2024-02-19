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

using namespace IMP;

int main() {

    try {

        Timer timer;

        std::cout << "***Tetrahedral to polyhedral conversion program***\n";

        Mesh<3,4> mesh;
        mesh.LoadFrom("/home/mood/Dropbox/phynetouch_sims/meshes/3D/inp/embc/abbott_electrode_filled.2mm.1R0d.inp");

        TetToPoly polygen;
        polygen.SetFeatureAngle(40.);
        std::cout << "***Primal mesh construction***\n\n";
        std::cout << "Loading primal mesh... ";
        polygen.LoadPrimalMesh(mesh);
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();

        std::cout << "Orient primal cells... ";

        polygen.OrientPrimalCells();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();
        std::cout << "Extract primal faces...";
        polygen.ExtractPrimalFaces();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();
        std::cout << "Extract primal edges...";
        polygen.ExtractPrimalEdges();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();
        std::cout << "Find edge features...";
        polygen.FindEdgesOnFeature();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();
        std::cout << "Find boundary verts...";
        polygen.FindVertsOnBoundaryAndFeature();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();

        std::cout << "\n***Dual mesh generation***\n\n";
        std::cout << "Generate dual verts...";
        polygen.GenerateDualVerts();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();
        std::cout << "Generate dual faces...";
        polygen.GenerateDualFaces();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();

        std::cout << "Generate dual cells...";
        polygen.GenerateDualCells();
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();

        std::cout << "Report...\n";
        std::cout << "cells: " << polygen.DualCells().size() << std::endl;
        std::cout << "facets: " << polygen.DualFaces().size() << std::endl;

        std::cout << "Save dual mesh data...";
        polygen.SaveDualCells("/home/mood/Dropbox/phynetouch_sims/meshes/3D/vtk/embc/abbott_electrode_filled.2mm.1R0d.evt.cells.vtk");
        polygen.SaveDualMesh( "/home/mood/Dropbox/phynetouch_sims/meshes/3D/evt/embc/abbott_electrode_filled.2mm.1R0d.evt");
        std::cout << timer.PrintElapsedTime() << std::endl;
        timer.Reset();

        std::cout << "Main program terminated sucessfully\n";

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
