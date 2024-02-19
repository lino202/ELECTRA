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

#ifndef IMP_ENGINE_TESSELATIONS_VORONOI_TPP_
#define IMP_ENGINE_TESSELATIONS_VORONOI_TPP_

#include "IMP/engine/tesselations/voronoi.hpp"


namespace IMP {

template <short DIM>
Voronoi<DIM>::Voronoi() : nodes_(), points_(), cells_(), facets_(), facets_on_nodes_(), cell_measures_(), node_sets_()
{}


template <short DIM>
Voronoi<DIM>::~Voronoi()
{}


template <short DIM>
void Voronoi<DIM>::LoadFrom(const std::string &filename)
{
    // Check if filename is given.
    if (filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not load voronoi tesselation file. The given filename is empty."));
    }

    // Get the extension of the mesh filename.
    auto ext = filename.substr(filename.length()-4);

    // Clear the voronoi tesselation containers.
    this->nodes_.clear();
    this->points_.clear();
    this->cells_.clear();
    this->facets_.clear();
    this->node_sets_.clear();

    // Load voronoi tesselation of given format.
    if (ext == ".evt") {
        EvtIO<DIM> evt_io;
        evt_io.ReadVoronoiFrom(filename);
        evt_io.LoadNodesIn(this->nodes_);
        evt_io.LoadPointsIn(this->points_);
        evt_io.LoadFacetsIn(this->facets_);
        evt_io.LoadCellsIn(this->cells_);
        if (evt_io.HasNodeSets()) {
            evt_io.LoadNodeSetsIn(this->node_sets_);
        }

        // Set parent and neighbor cell ids to facets.
        int cid = 0;
        for (const auto &cell : this->cells_) {
            for (const auto &fid : cell.Connectivity()) {
                (this->facets_[fid].ParentCellId() == -1) ?
                    this->facets_[fid].SetParentCellId(cid) : this->facets_[fid].SetNeighCellId(cid);
            }
            cid++;
        }

        // Set for each node its children facets.
        this->facets_on_nodes_.clear();
        this->facets_on_nodes_.resize(this->nodes_.size());
        for (std::size_t fid = 0; fid != this->facets_.size(); ++fid) {
            this->facets_on_nodes_[this->facets_[fid].ParentCellId()].emplace_back(fid);
        }

    } else {
        std::string err_msg = "Could not load voronoi tesselation of unknown format. Check given filename: " + filename;
        throw std::invalid_argument(Logger::Error(err_msg));
    }
}


template <short DIM>
void Voronoi<DIM>::SaveTo(const std::string &filename)
{
    // Check if mesh filename is given.
    if (filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not save voronoi tesselation to file. The given filename is empty."));
    }

    // Get the extension of the mesh filename.
    auto ext = filename.substr(filename.length()-4);

    // Load voronoi tesselation of given format.
    if (ext == ".evt") {
        EvtIO<DIM> evt_io;
        evt_io.SaveVoronoiTo(filename, this->nodes_, this->node_sets_, this->points_, this->facets_, this->cells_);
    } else {
        std::string err_msg = "Could not save voronoi tesselation to unknown format. Check given filename: " + filename;
        throw std::invalid_argument(Logger::Error(err_msg));
    }
}


template <short DIM>
void Voronoi<DIM>::Scale(double scale_value)
{
    for (auto &node : this->nodes_) {
        node *= scale_value;
    }
}


template <short DIM>
void Voronoi<DIM>::AddNodeSet(const NodeSet &nodeset)
{
    this->node_sets_[nodeset.Name()] = nodeset;
}


template <short DIM>
void Voronoi<DIM>::ComputeCellMeasures()
{
    if (this->cells_.size() == 0 || this->points_.size() == 0 || this->facets_.size() == 0) {
        auto err_msg = "Could not compute voronoi cell measures. Set cells, points, and facets first.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    this->cell_measures_.clear();
    this->cell_measures_.resize(this->cells_.size());
    for (auto i = 0; i != this->CellsNum(); ++i)
        this->cell_measures_[i] = this->cells_[i].Measure(this->points_, this->facets_);
}


template <short DIM>
Vec<DIM,double> Voronoi<DIM>::CellCentroid(std::size_t id) const
{
    return this->cells_[id].Centroid(this->points_, this->facets_);
}


} // End of namespace IMP.

#endif // IMP_ENGINE_TESSELATIONS_VORONOI_TPP_
