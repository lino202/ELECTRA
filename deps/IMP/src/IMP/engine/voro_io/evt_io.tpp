/*
 * IMP. Image and Mesh Processing library.
 * Copyright (C) 2016  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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

#ifndef IMP_ENGINE_VORO_IO_EVT_IO_TPP_
#define IMP_ENGINE_VORO_IO_EVT_IO_TPP_

#include "IMP/engine/voro_io/evt_io.hpp"

namespace IMP {

template <short DIM>
EvtIO<DIM>::EvtIO() : parsed_voronoi_(), mapped_node_ids_(), mapped_point_ids_(), mapped_facet_ids_(), nsets_startlines_(),
    nodes_startline_(0), points_startline_(0), facets_startline_(0), cells_startline_(0), has_nsets_(false)
{}


template <short DIM>
EvtIO<DIM>::~EvtIO()
{}


template <short DIM>
void EvtIO<DIM>::ReadVoronoiFrom(const std::string &filename)
{
    // Clear containers.
    this->parsed_voronoi_.clear();
    this->nsets_startlines_.clear();

    // Check if file is in ELECTRA Voronoi Tesselation format (.evt).
    if (fs::path(filename).extension().string() != ".evt") {
        auto err_msg = "Expected file in ELECTRA Voronoi Tesselation format (*.evt)";
        throw std::invalid_argument(Logger::Error(err_msg));
    }

    // Open file.
    auto file = std::ifstream(filename, std::ios::in);
    if (!file.is_open()) {
        auto err_msg = "Could not open the voronoi tesselation file.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Available voronoi information conditionals to default false values.
    bool  has_nodes   = false;
    bool  has_points  = false;
    bool  has_cells   = false;
    bool  has_facets  = false;
    this->has_nsets_  = false;

    //Read voronoi file line by line.
    auto line = std::string{""};
    auto line_number = int{0};
    while (std::getline(file, line)) {

        // Read mesh file lines.
        this->parsed_voronoi_.emplace_back(line);

        // Check if nodes list was found in the mesh.
        if (ALGORITHMS::ExistWord(line, "*nodes")) {
            this->nodes_startline_ = line_number;
            has_nodes = true;
        }

        // Check if points list was found in the mesh.
        if (ALGORITHMS::ExistWord(line, "*points")) {
            this->points_startline_ = line_number;
            has_points = true;
        }

        // Check if facets list was found in the mesh.
        if (ALGORITHMS::ExistWord(line, "*facets")) {
            this->facets_startline_ = line_number;
            has_facets = true;
        }

        // Check if cells list was found in the mesh.
        if (ALGORITHMS::ExistWord(line, "*cells")) {
            this->cells_startline_ = line_number;
            has_cells = true;
        }

        // Find the first lines of the boundary node sets.
        if (ALGORITHMS::ExistWord(line, "*nset")) {
            this->nsets_startlines_.emplace_back(line_number);
            if (!this->has_nsets_) { this->has_nsets_ = true; }
        }

        // Increase line number.
        line_number++;
    }

    // Close the voronoi file.
    file.close();

    // Check if mandatory information is available in the voronoi tesselation file.
    if (!has_nodes && !has_points && !has_facets && !has_cells) {
        auto err_msg = "Something went wrong while reading nodes, points, facets, or cells. Check voronoi tesselation file: " + filename;
        throw std::runtime_error(Logger::Error(err_msg));
    }

}


template <short DIM>
void EvtIO<DIM>::SaveVoronoiTo(const std::string &filename, const std::vector<Vec<DIM, double>> &nodes,
    const std::unordered_map<std::string, NodeSet> &node_sets, const std::vector<Vec<DIM, double>> &points,
    const std::vector<PolyFacet> &facets, const std::vector<PolyCell> &cells) const
{
    // Check if filename is given.
    if (filename.empty()) {
        auto err_msg = "Could not save the voronoi tesselation. Check given filename: " + filename;
        throw std::invalid_argument(Logger::Error(err_msg));
    }

    // Open output file.
    auto out = std::ofstream{};
    if (fs::path(filename).extension().string() != ".evt") {
        out.open(filename+".evt", std::ios::out);
    } else {
        out.open(filename, std::ios::out);
    }
    if (!out.is_open()) {
        auto err_msg = "Could not create the voronoi tesselation output file. Check given filename: " + filename;
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Save header information.
    out << "*\n";
    out << "** ELECTRA Voronoi Tesselation file\n";
    out << "** Generated by: IMP library\n";
    if (DIM==3) out << "** Cell connectivity corresponds to the connected facets that compose the cell\n";
    else out << "** Cell connectivity corresponds to the connected points that compose the cell\n";
    out << "**\n";

    // Save voronoi tesselation nodes coordinates.
    out << "*Nodes " << nodes.size() << "\n";
    for (const auto &node : nodes) {
        // Save node id.
        out << &node - &nodes[0] + 1 << ", ";

        // Save node coordinates.
        for (short i = 0; i != DIM-1; ++i) {
            out << std::setprecision(15) << node[i] << ", ";
        }
        out << node[DIM-1] << "\n";
    }

    // Save voronoi tesselation points coordinates.
    out << "*Points " << points.size() << "\n";
    auto pid = int{0};
    for (const auto &point : points) {
        // Save point id.
        out << ++pid << ", ";

        // Save point coordinates.
        for (short i = 0; i != DIM-1; ++i) {
            out << std::setprecision(15) << point[i] << ", ";
        }
        out << point[DIM-1] << "\n";
    }

    // Save voronoi tesselation facets if available.
    out << "*Facets " << facets.size() << "\n";
    auto fid = int{0};
    for (const auto &facet : facets) {
        // Save facet id.
        out << ++fid << ", " << facet.Connectivity().size() << ", ";

        // Save facet connectivity.
        for (std::size_t i = 0; i != facet.Connectivity().size()-1; ++i) out << facet.C(i)+1 << ", ";
        out << facet.C(facet.Connectivity().size()-1)+1 << "\n";
    }

    // Save voronoi tesselation cells.
    out << "*Cells " << cells.size() << "\n";
    auto cid = int{0};
    for (const auto &cell : cells) {

        // Save cell id.
        out << ++cid << ", " << cell.Connectivity().size() << ", ";

        // Save cell connectivity.
        for (std::size_t i = 0; i != cell.Connectivity().size()-1; ++i) out << cell.C(i)+1 << ", ";
        out << cell.C(cell.Connectivity().size()-1)+1 << "\n";
    }

    // Save voronoi tesselation nodesets.
    for (const auto &nset : node_sets) {
        out << "*Nset, nset=" << nset.first << "\n";
        auto count = std::size_t{0}, pos = std::size_t{0};
        for (const auto  &nid : nset.second.NodeIds()) {
            pos = &nid - &nset.second.NodeIds()[0];
            out << nid+1;
            count++;
            if (count < 15 && pos != nset.second.NodeIds().size()-1) { out << ", "; }
            else { out << "\n"; count = 0; }
        }
    }

    // Save ending header.
    out << "*End\n" << "**\n";

    // Close the output file.
    out.close();

}


template <short DIM>
void EvtIO<DIM>::LoadNodesIn(std::vector<Vec<DIM, double>> &nodes)
{
    // Clean the mapped node ids container.
    this->mapped_node_ids_.clear();

    // Clean the nodes container and reserve space.
    int nodes_num = 0;
    std::string line = this->parsed_voronoi_[this->nodes_startline_];
    std::stringstream ss(line);
    ss >> line >> nodes_num;
    ss.str(std::string{});  ss.clear();
    nodes.clear();
    nodes.resize(nodes_num, Vec<DIM,double>{});

    // The coordinates of a node in the voronoi tesselation.
    Vec<DIM, double> coords;

    // Iterate though voronoi tesselation starting from the nodes list starting line.
    // Skip the first line to start from the first node's coordinates.
    std::size_t key_id = 0, cnt = 0;
    for (std::size_t it = this->nodes_startline_+1; it != this->parsed_voronoi_.size(); ++it) {

        // Get line and replace comma occurences with space for stringstream processing.
        line = this->parsed_voronoi_[it];
        std::replace(line.begin(), line.end(), ',', ' ');

        // Get the coordinates until reach the end of the nodes list.
        ss.str(line);
        if (!(ss >> key_id >> coords)) { break; }
        ss.str(std::string{});
        ss.clear();

        // Store the current node's coordinates.
        nodes[cnt++] = coords;

        // Store the mapping from the key id of the current node to the container's id.
        this->mapped_node_ids_[key_id] = it - (this->nodes_startline_+1);
    }

}

template <short DIM>
void EvtIO<DIM>::LoadPointsIn(std::vector<Vec<DIM, double>> &points)
{

    // Clean the mapped node ids container.
    this->mapped_point_ids_.clear();

    // Clean the points container and reserve space.
    int points_num = 0;
    std::string line = this->parsed_voronoi_[this->points_startline_];
    std::stringstream ss(line);
    ss >> line >> points_num;
    ss.str(std::string{}); ss.clear();
    points.clear();
    points.resize(points_num, Vec<DIM,double>{});

    // The coordinates of a point in the voronoi tesselation.
    Vec<DIM, double> coords;

    // Iterate though voronoi tesselation starting from the points list starting line.
    // Skip the first line to start from the first point's coordinates.
    std::size_t key_id = 0, cnt = 0;
    for (std::size_t it = this->points_startline_+1; it != this->parsed_voronoi_.size(); ++it) {

        // Get line and replace comma occurences with space for stringstream processing.
        line = this->parsed_voronoi_[it];
        std::replace(line.begin(), line.end(), ',', ' ');

        // Get the coordinates until reach the end of the points list.
        ss.str(line);
        if (!(ss >> key_id >> coords)) { break; }
        ss.str(std::string{});
        ss.clear();

        // Store the current point's coordinates.
        points[cnt++] = coords;

        // Store the mapping from the key id of the current point to the container's id.
        this->mapped_point_ids_[key_id] = it - (this->points_startline_+1);
    }
}


template <short DIM>
void EvtIO<DIM>::LoadFacetsIn(std::vector<PolyFacet> &facets)
{
    // Check if mapped point ids have been set.
    if (this->mapped_point_ids_.empty()) {
        throw std::runtime_error(Logger::Error("Could not load voronoi tesselation facets. Load the voronoi tesselation points first."));
    }

    // Get number of facets.
    int facets_num = 0;
    std::string line = this->parsed_voronoi_[this->facets_startline_];
    std::stringstream ss(line);
    ss >> line >> facets_num;
    ss.str(std::string{});  ss.clear();

    // Clean the facets container.
    facets.clear();
    facets.resize(facets_num, PolyFacet{});

    // Clean the mapped facet ids container.
    this->mapped_facet_ids_.clear();
    this->mapped_facet_ids_.reserve(facets_num);

    // Iterate though the voronoi tesselation starting from the first line of the facets list.
    // Skip a line to read the first facets's connectivity.
    PolyFacet facet;
    int key_id = 0, facet_size = 0, conn_id = 0, cnt = 0;
    for (std::size_t it = this->facets_startline_+1; it != this->parsed_voronoi_.size(); ++it) {

        // Get line and replace comma occurences with space for stringstream processing.
        line = this->parsed_voronoi_[it];
        std::replace(line.begin(), line.end(), ',', ' ');

        // Get the key id and the size of the current facet.
        ss.str(line);
        if (!(ss >> key_id >> facet_size)) { break; }

        // Assign connectivity to the corresponding facet in the container.
        // Correct for any offset in points container.
        int i = 0;
        while (i != facet_size) {
            ss >> conn_id;
            facets[cnt].AddInConnectivity(this->mapped_point_ids_.at(conn_id));
            i++;
        }
        cnt++;
        ss.str(std::string{});
        ss.clear();

        // Store the mapping from the key id of the current cell to the container's id.
        this->mapped_facet_ids_[key_id] = it - (this->facets_startline_+1);
    }
}


template <short DIM>
void EvtIO<DIM>::LoadCellsIn(std::vector<PolyCell> &cells)
{
    // Check if mapped facets or point ids have been set.
    if (this->mapped_facet_ids_.empty()) {
        throw std::runtime_error(Logger::Error("Could not load voronoi tesselation cells. Load the voronoi tesselation facets first."));
    }

    // Clean the cells container.
    cells.clear();

    // Clean the cells container and reserve space.
    int cells_num = 0;
    std::string line = this->parsed_voronoi_[this->cells_startline_];
    std::stringstream ss(line);
    ss >> line >> cells_num;
    ss.str(std::string{});  ss.clear();
    cells.clear();
    cells.resize(cells_num, PolyCell{});

    // Iterate though the voronoi tesselation starting from the first line of the cells list.
    // Skip a line to read the first cell's connectivity.
    PolyCell cell;
    int key_id = 0, cell_size = 0, conn_id = 0, cnt = 0;
    for (std::size_t it = this->cells_startline_+1; it != this->parsed_voronoi_.size(); ++it) {

        // Get line and replace comma occurences with space for stringstream processing.
        line = this->parsed_voronoi_[it];
        std::replace(line.begin(), line.end(), ',', ' ');

        // Get the key id and the size of the current cell.
        ss.str(line);
        if (!(ss >> key_id >> cell_size)) { break; }

        // Assign connectivity to the corresponding cell in the container.
        // Correct for any offset in points container.
        int i = 0;
        while (i != cell_size) {
            ss >> conn_id;
            cells[cnt].AddInConnectivity(this->mapped_facet_ids_.at(conn_id));
            i++;
        }
        cnt++;
        ss.str(std::string{});
        ss.clear();
    }
}


template <short DIM>
void EvtIO<DIM>::LoadNodeSetsIn(std::unordered_map<std::string, NodeSet> &node_sets)
{
    // Initialize the node_sets container.
    node_sets.clear();

    // The node's index. Initialized in invalid value (-1).
    int node_id = -1;

    // Create a temporary node set.
    NodeSet nset;

    // Iterate though mesh starting from each node set's starting line.
    for (auto &current_nset_startline: this->nsets_startlines_) {

        // Clear temporary node set.
        nset.Clear();

        // Extract node set name from its starting line.
        std::string nset_name = this->parsed_voronoi_.at(current_nset_startline);
        nset_name = nset_name.substr(nset_name.find_last_of("=")+1);
        nset_name.erase(std::remove_if(nset_name.begin(), nset_name.end(), isspace), nset_name.end());

        // Assign the node set name to the temporary node set.
        nset.SetName(nset_name);

        // Iterate node indices in the node set.
        // Skipping the node set title line.
        std::string line = "";
        std::stringstream ss;
        for (std::vector<std::string>::size_type it = current_nset_startline+1;
                                                 it != this->parsed_voronoi_.size(); ++it) {

            // Get line and replace comma occurences with space for stringstream processing.
            line = this->parsed_voronoi_[it];
            std::replace(line.begin(), line.end(), ',', ' ');

            // Convert line to stringstream to be passed in the connectivity variables.
            ss.str(line);

            // Stop loop if header line found.
            if (line.find('*') != std::string::npos) { break; }

            // Store the indices of the nodes belonging to the nodeset.
            while (ss >> node_id) {
                // Add the mapped index of the node id to the temporary node set.
                nset.EditNodeIds().emplace_back(this->mapped_node_ids_[node_id]);
            }

            ss.str(std::string{});
            ss.clear();

        } // End of Iterate node indices in the node set.

        // Assing temporary node set in the node sets container.
        node_sets[nset_name] = nset;

    } // End of Iterate though mesh starting from each node set's starting line.

}


}  //end of namespace IMP

#endif //IMP_ENGINE_VORO_IO_EVT_IO_TPP_
