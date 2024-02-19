/*
 * IMP. Image to Mesh Processing library.
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


/**
   \file grid.hpp
   \brief Grid class header file.
   \author Konstantinos A. Mountris
   \date 21/05/2019
*/

#ifndef IMP_ENGINE_TESSELATIONS_GRID_HPP_
#define IMP_ENGINE_TESSELATIONS_GRID_HPP_

#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/elements/cell.hpp"
#include "IMP/engine/elements/cell_props.hpp"
#include "IMP/engine/mesh_io/abaqus_io.hpp"
#include "IMP/engine/mesh_io/vtk_io.hpp"
#include "IMP/engine/topology/node_set.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <Eigen/Sparse>

#include <vector>
#include <unordered_map>

namespace IMP {

/** \addtogroup Tesselations \{ */

/**
 * \class Grid
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a grid of nodes without any connectivity information.
 * \tparam DIM the spatial dimensions of the grid.
 */
template<short DIM, short CELL_VERTS=1>
class Grid {

private:

    std::vector<Vec<DIM, double>> nodes_;                           /**< The nodes of the grid */

    std::vector<Cell<DIM, CELL_VERTS>>  ghost_cells_;               /**< The ghost cells of the grid */

    std::unordered_map<std::string, NodeSet> node_sets_;            /**< The node sets of the grid */

    Eigen::SparseMatrix<double, Eigen::RowMajor> node_normals_;     /**< The normals of the grid nodes */

    CellShape cells_shape_;                                         /**< The shape of the ghost cells of the grid */


protected:

    /**
     * \brief Update the type of the ghost cells of the grid.
     * \return [void] 
     */
    inline void UpdateGhostCellsShape();


public:
    /**
     * \brief Grid constructor.
     */
    inline Grid();


    /**
     * \brief Grid destructor.
     */
    inline virtual ~Grid();


    /**
     * \brief Load the grid from a file.
     * \param [in] grid_filename The name of the grid's file.
     * \return [void]
     */
    inline void LoadFrom(const std::string &grid_filename);


    /**
     * \brief Load the normal vectors for the surfaces grid nodes from a file.
     * \param [in] normals_filename The name of the file with the normal vectors.
     * \return [void] 
     */
    inline void LoadNodeNormals(const std::string &normals_filename);


    /**
     * \brief Write access to the nodes of the grid.
     * \return [std::vector<IMP::Vec<DIM, double> >&] The nodes of the grid with write access.
     */
    inline std::vector<Vec<DIM, double>> & EditNodes() { return this->nodes_; }


    /**
     * \brief Read-only access to the nodes of the grid.
     * \return [const std::vector<IMP::Vec<DIM, double> >&] The nodes of the grid nodes with read-only access.
     */
    inline const std::vector<Vec<DIM, double>> & Nodes() const { return this->nodes_; }


    /**
     * \brief Get a single node of the grid.
     * Fast access, no range check. 
     * \param [in] id The index of the node. 
     * \return [const Vec<DIM, double>&] The node with index id. 
     */
    inline const Vec<DIM, double> & Nodes(std::size_t id) const { return this->nodes_[id]; }


    /**
     * \brief Get a single node of the grid.
     * Slow access, with range check.
     * \param [in] id The index of the node. 
     * \return [const Vec<DIM, double>&] The node with index id. 
     */
    inline const Vec<DIM, double> & NodesAt(std::size_t id) const { return this->nodes_.at(id); }


    /**
     * \brief Get the number of nodes in the grid.
     * \return [int] The number of nodes in the grid. 
     */
    inline int NodesNum() const { return static_cast<int>(this->nodes_.size()); }


    /**
     * \brief Write access to the ghost cells of the grid.
     * \return [std::vector<Cell<DIM, CELL_VERTS>>&] The ghost cells of the grid with write access.
     */
    inline std::vector<Cell<DIM, CELL_VERTS>> & EditGhostCells() { return this->ghost_cells_; }


    /**
     * \brief Read-only access to the ghost cells of the grid.
     * \return [const std::vector<Cell<DIM, CELL_VERTS>>&] The ghost cells of the grid with read-only access.
     */
    inline const std::vector<Cell<DIM, CELL_VERTS>> & GhostCells() const { return this->ghost_cells_; }


    /**
     * \brief Get a single ghost cell of the grid.
     * 
     * Fast access, no range check.
     * 
     * \param [in] id The index of the ghost cell. 
     * \return [const Cell<DIM, CELL_VERTS>&] The ghost cell with index id. 
     */
    inline const Cell<DIM, CELL_VERTS> & GhostCells(std::size_t id) const { return this->ghost_cells_[id]; }


    /**
     * \brief Get a single ghost cell of the grid.
     * 
     * Slow access, with range check.
     * 
     * \param [in] id The index of the ghost cell. 
     * \return [const Vec<DIM, double>&] The ghost cell with index id. 
     */
    inline const Cell<DIM, CELL_VERTS> & GhostCellsAt(std::size_t id) const { return this->ghost_cells_.at(id); }


    /**
     * \brief Get the number of the ghost cells in the grid.
     * \return [int] The number of the ghost cells in the grid. 
     */
    inline int GhostCellsNum() const { return static_cast<int>(this->ghost_cells_.size()); }


    /**
     * \brief Get the node sets of the grid.
     * The node sets are groups of nodes sharing common attributes.
     * \return [const std::unordered_map<std::string, IMP::NodeSet>&] The node sets of the grid.
     */
    inline const std::unordered_map<std::string, NodeSet> & NodeSets() const { return this->node_sets_; }


    /**
     * \brief Get the node set of the grid that is specified by name.
     * \param [in] name The name of the node set. 
     * \return [const IMP::NodeSet&] The node set of the grid that is specified by name. 
     */
    inline const NodeSet & NodeSets(const std::string &name) const { return this->node_sets_.at(name); }


    /**
     * \brief Get the number of node sets of the grid. 
     * \return [int] The number of node sets of the grid.
     */
    inline int NodeSetsNum() const { return static_cast<int>(this->nodes_sets_.size()); }


    /**
     * \brief Get the normals of the grid nodes.
     * \return [const Eigen::SparseMatrix<double, Eigen::RowMajor>&] The normals of the grid nodes.
     */
    inline const Eigen::SparseMatrix<double, Eigen::RowMajor> & NodeNormals() const { return this->node_normals_; }

};


/** \} End of Doxygen Groups*/

} // End of namespace IMP

#endif //IMP_ENGINE_TESSELATIONS_GRID_HPP_

#include "IMP/engine/tesselations/grid.tpp"