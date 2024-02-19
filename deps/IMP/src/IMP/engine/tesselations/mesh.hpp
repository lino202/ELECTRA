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
   \file mesh.hpp
   \brief Mesh class header file.
   \author Konstantinos A. Mountris
   \date 11/10/2018
*/

#ifndef IMP_ENGINE_TESSELATIONS_MESH_HPP_
#define IMP_ENGINE_TESSELATIONS_MESH_HPP_

#include "IMP/engine/elements/cell.hpp"
#include "IMP/engine/elements/cell_props.hpp"
#include "IMP/engine/elements/half_facet.hpp"
#include "IMP/engine/elements/half_facets_array.hpp"
#include "IMP/engine/mesh_io/abaqus_io.hpp"
#include "IMP/engine/mesh_io/vtk_io.hpp"
#include "IMP/engine/topology/node_set.hpp"
#include "IMP/engine/utilities/logger.hpp"
#include "IMP/engine/utilities/timer.hpp"

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <utility>

namespace IMP {

/** \addtogroup Tesselations \{ */

/**
 * \class Mesh
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a generic mesh.
 * 
 * The mesh can be of any dimension and composed of any elements of the same type. 
 * It is based on the Array-based Half-Facets structure for efficient topological querries.
 * 
 * \tparam DIM The geometrical dimension of the mesh.
 * \tparam CELL_VERTS_NUM The number of vertices per cell.
 */
template <short DIM, short CELL_VERTS>
class Mesh {

private:

    std::vector<Vec<DIM, double> > nodes_;                  /**< The nodes of the mesh */

    std::vector<Cell<DIM, CELL_VERTS> >  cells_;            /**< The cells of the mesh */

    std::unordered_map<std::string, NodeSet> node_sets_;    /**< The node sets of the mesh */

    HalfFacetsArray nodes_to_half_facets_;                  /**< The mapping of the mesh nodes to an incident half-facet */

    CellShape cells_shape_;                                 /**< The shape of the cells in the mesh */

    bool are_sibling_half_facets_set_;                      /**< Boolean to check if sibling half facets have been set in the mesh cells */

    bool are_nodes_mapped_to_half_facets_;                  /**< Boolean to check if nodes have been mapped to the half facets */


protected:

    /**
     * \brief Update the shape of the cells in the mesh.
     * Cells with different shapes are not supported. 
     * \return [void]
     */
    inline void UpdateCellsShape();


public:

    /**
     * \brief Mesh constructor.
     */
    inline Mesh();


    /**
     * \brief Mesh destructor.
     */
    inline virtual ~Mesh();


    /**
     * \brief Load a mesh from a file.
     * The following mesh file formats are supported: *.inp
     * \param [in] mesh_filename The filename of the mesh file to be loaded.
     */
    inline void LoadFrom(const std::string &mesh_filename, bool built_ahf = true);


    /**
     * \brief Save a mesh to a file.
     * The following mesh file formats are supported: *.inp
     * \param [in] mesh_filename The filename of the mesh file to be saved.
     */
    inline void SaveTo(const std::string &mesh_filename);


    /**
     * @brief 
     * @param scale_value 
     */
    inline void Scale(double scale_value);



    /**
     * \brief Append a new node in the mesh. 
     * 
     * \param [in] coords The node's coordinates.
     * \return [void] 
     */
    inline void AppendNode(std::initializer_list<double> coords);


    /**
     * \brief Append a new cell in the mesh. 
     * 
     * \param connectivity The nodes connectivity of the cell. 
     */
    inline void AppendCell(std::initializer_list<int> connectivity);


    /**
     * \brief Add a new nodeset in the mesh.
     * \param [in] nodeset The new nodeset to be added.
     * \return [void]
     */
    inline void AddNodeSet(const NodeSet &nodeset);


    /**
     * \brief
     */
    inline void BuildAhf() { this->FindSiblingHalfFacets(); this->MapNodesToHalfFacets(); }


    /**
     * \brief Find the cells that are attached to a given node.
     * 
     * \param [in] node_id 
     * \return std::vector<int> 
     */
    inline std::vector<int> AttachedCellIdsToNode(int node_id) const;


    /**
     * \brief Find the cells that are attached to a given node and additionally share a half-facet with the given cell.
     * 
     * \param [in] node_id 
     * \param [in] cell_id 
     * \return std::vector<int> 
     */
    inline std::vector<int> AttachedCellIdsToNodeAndCell(int node_id, int cell_id) const;


    /**
     * \brief Get the half facets of free boundary of the mesh.
     * \return [std::vector<HalfFacet>] The half facets on the free boundary of the mesh.
     */
    inline std::vector<HalfFacet> FreeFacets() const;


    /**
     * \brief Get the indices of the nodes on the free boundary of the mesh.
     * \return [std::vector<int>] The indices of the nodes on the free boundary of the mesh.
     */
    inline std::vector<int> FreeNodeIds() const;


    /**
     * \brief Get the nodes of the mesh.
     * \return [const std::vector<IMP::Vec<DIM, double> >&] The nodes of the mesh.
     */
    inline const std::vector<Vec<DIM, double>> & Nodes() const { return this->nodes_; }


    /**
     * \brief Get a single node of the mesh.
     * Fast access, no range check.
     * \param [in] id The index of the node. 
     * \return [const Vec<DIM, double>&] The node with index id. 
     */
    inline const Vec<DIM, double> & Nodes(std::size_t id) const { return this->nodes_[id]; }


    /**
     * \brief Get a single node of the mesh.
     * Slow access, with range check.
     * \param [in] id The index of the node. 
     * \return [const Vec<DIM, double>&] The node with index id. 
     */
    inline const Vec<DIM, double> & NodesAt(std::size_t id) const { return this->nodes_.at(id); }


    /**
     * \brief Get the number of nodes in the mesh.
     * \return [int] The number of nodes in the mesh. 
     */
    inline int NodesNum() const { return static_cast<int>(this->nodes_.size()); }


    /**
     * \brief Get the cells of the mesh.
     * \return [const std::vector<IMP::Cell<DIM, CELL_VERTS> >&] The cells of the mesh.
     */
    inline const std::vector<Cell<DIM, CELL_VERTS>> & Cells() const { return this->cells_; }


    /**
     * \brief Get a single cell of the mesh.
     * Fast access, no range check.
     * \param [in] id The index of the cell. 
     * \return [const Cell<DIM, CELL_VERTS>&] The cell with index id. 
     */
    inline const Cell<DIM, CELL_VERTS> & Cells(std::size_t id) const { return this->cells_[id]; }


    /**
     * \brief Get a single cell of the mesh.
     * Slow access, with range check.
     * \param [in] id The index of the cell. 
     * \return [const Vec<DIM, double>&] The cell with index id. 
     */
    inline const Cell<DIM, CELL_VERTS> & CellsAt(std::size_t id) const { return this->cells_.at(id); }


    /**
     * \brief Get the number of cells in the mesh. 
     * \return [int] The number of cells in the mesh.
     */
    inline int CellsNum() const { return static_cast<int>(this->cells_.size()); }


    /**
     * \brief Get the type of the mesh cells shape.
     * \return [IMP::CellShape] The type of the mesh cells shape.
     */
    inline CellShape CellsShape() const { return this->cells_shape_; }


    /**
     * \brief Get the node sets of the mesh.
     * The node sets are groups of nodes sharing common attributes.
     * \return [const std::unordered_map<std::string, IMP::NodeSet>&] The node sets of the mesh.
     */
    inline const std::unordered_map<std::string, NodeSet> & NodeSets() const { return this->node_sets_; }


    /**
     * \brief Get the node set of the mesh that is specified by name.
     * \param [in] name The name of the required node set. 
     * \return [const NodeSet&] The node set of the mesh that is specified by name. 
     */
    inline const NodeSet & NodeSets(const std::string &name) const { return this->node_sets_.at(name); }


    /**
     * \brief Get the number of node sets of the mesh. 
     * \return [int] The number of node sets of the mesh.
     */
    inline int NodeSetsNum() const { return static_cast<int>(this->node_sets_.size()); }


    /**
     * \brief 
     * \return const HalfFacetsArray& 
     */
    inline const HalfFacetsArray & NodesToHalfFacets() const { return this->nodes_to_half_facets_; }


    /**
     * \brief Check if the AHF (Array-based HalfFacet) structure is built.
     * 
     * The AHF structure is considered built when both the sibling half facets of the mesh cells have been extracted
     * and the mesh nodes have been mapped to the half facets.
     * 
     * \return [true] The AHF structure is built.
     * \return [false] The AHF structure is not built.
     */
    inline bool IsAhfBuilt() const { return this->are_sibling_half_facets_set_ && this->are_nodes_mapped_to_half_facets_; }


// protected:

    /**
     * \brief Find and assign the sibling half-facets to each cell of the mesh.
     * \return [void] 
     */
    inline void FindSiblingHalfFacets();


    /**
     * \brief Map the mesh nodes to an incident half-facet. Boundary half-facets are prioritized.
     * \return [void]
     */
    inline void MapNodesToHalfFacets();


};



/** \} End of Doxygen Groups*/
} //namespace IMP

#include "IMP/engine/tesselations/mesh.tpp"

#endif //IMP_ENGINE_TESSELATIONS_MESH_HPP_
