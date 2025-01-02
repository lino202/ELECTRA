/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019
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

#ifndef ELECTRA_APPS_TOOLS_CONFIG_GEO_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_GEO_TPP_

#include "config_geo.hpp"

namespace APP_ELECTRA
{

template<short DIM, short CELL_NODES>
ConfigGeo<DIM, CELL_NODES>::ConfigGeo()
{}


template<short DIM, short CELL_NODES>
ConfigGeo<DIM, CELL_NODES>::~ConfigGeo()
{}


template<short DIM, short CELL_NODES>
void ConfigGeo<DIM, CELL_NODES>::SetGeometryModel(const Parser &parser, IMP::Mesh<DIM, CELL_NODES> &mesh, IMP::Grid<DIM, CELL_NODES> &grid,
                                 IMP::Voronoi<DIM> &voro, std::ostream &stream) const
{
    // Resolve geometry file path.
    std::string filename = parser.GetValue<std::string>("tissue.geometry.file");
    std::filesystem::path filepath(filename);
    std::string parent_path_string{filepath.parent_path()};
    if (!filepath.has_root_path() || parent_path_string.size() == 1) {
        if (filename[0] == '/' || filename[0] == '\\') {
            filename = parser.ParentPath() + filename;
        } else {
            filename = parser.ParentPath() + "/" + filename;
        }
    }

    // Get type of numerical method.
    std::string method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);

    // Set geometry according to the numerical method.
    if (method == "fem") {
        mesh.LoadFrom(filename);

        stream << ELECTRA::Logger::Message("Loaded mesh: " + filename + "\n");
        stream << ELECTRA::Logger::Message("Number of nodes: ") << mesh.NodesNum() << "\n";
        stream << ELECTRA::Logger::Message("Number of cells: ") << mesh.CellsNum() << "\n";
    } else if (method == "fpm") {
        // Temporary load of mesh for output.
        mesh.LoadFrom(filename);

        std::string voro_file = parser.GetValue<std::string>("numerical approximation.fpm.voronoi");
        std::filesystem::path voro_fp(voro_file);
        std::string voro_fp_path_string{voro_fp.parent_path()};
        if (!voro_fp.has_root_path() || voro_fp_path_string.size() == 1) {
            if (voro_file[0] == '/' || voro_file[0] == '\\') {
                voro_file = parser.ParentPath() + voro_file;
            } else {
                voro_file = parser.ParentPath() + "/" + voro_file;
            }
        }
        voro.LoadFrom(voro_file);

        stream << ELECTRA::Logger::Message("Loaded mesh: " + filename + "\n");
        stream << ELECTRA::Logger::Message("Loaded Voronoi tesselation: " + voro_file + "\n");
        stream << ELECTRA::Logger::Message("Number of nodes: ") << voro.NodesNum() << "\n";
        stream << ELECTRA::Logger::Message("Number of Voronoi points: ") << voro.PointsNum() << "\n";
        if (voro.FacetsNum() != 0)
            stream << ELECTRA::Logger::Message("Number of Voronoi facets: ") << voro.FacetsNum() << "\n";
        stream << ELECTRA::Logger::Message("Number of Voronoi cells: ") << voro.CellsNum() << "\n";

    } else if (method == "mcm") {
        grid.LoadFrom(filename);

        // Check for surface nodes' normals file.
        if (!parser.HasAttribute("tissue.geometry.normals")) {
            std::string error_str = "tissue.geometry.normals attribute required for Mixed Collocation method simulation but not found.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }

        // Resolve surface nodes normals file path.
        std::string normals_filename = parser.GetValue<std::string>("tissue.geometry.normals");
        std::filesystem::path pn(normals_filename);
        std::string parent_path_string{pn.parent_path()};
        if (!pn.has_root_path() || parent_path_string.size() == 1) {
            if (normals_filename[0] == '/' || normals_filename[0] == '\\') {
                normals_filename = parser.ParentPath() + normals_filename;
            } else {
                normals_filename = parser.ParentPath() + "/" + normals_filename;
            }
        }

        // Load nodal normals.
        grid.LoadNodeNormals(normals_filename);
        stream << ELECTRA::Logger::Message("Loaded grid: " + filename + "\n");
        stream << ELECTRA::Logger::Message("Loaded surface node normals from file: " + normals_filename + "\n");
        stream << ELECTRA::Logger::Message("Number of nodes: ") << grid.NodesNum() << "\n";
    }
}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_UNITS_TPP_