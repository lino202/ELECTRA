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


/**
   \file vtk_io.hpp
   \brief VTK Simple Legacy format input/output class header file.
   \author Konstantinos A. Mountris
   \date 17/10/2019
*/

#ifndef IMP_ENGINE_MESH_IO_VTK_IO_HPP_
#define IMP_ENGINE_MESH_IO_VTK_IO_HPP_

#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/elements/cell.hpp"
#include "IMP/engine/topology/node_set.hpp"
#include "IMP/engine/elements/cell_props.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <iomanip>
#include <vector>
#include <unordered_map>
#include <map>
#include <iterator>
#include <utility>
#include <algorithm>
#include <string>

#include <exception>
#include <stdexcept>

#include <sstream>
#include <iostream>
#include <fstream>

namespace IMP {

/** \addtogroup MeshIO \{ */


/**
 * \brief Class implementing input/output for mesh in VTK Simple Legacy (*.vtk) format. 
 * 
 * \tparam DIM The dimensions of the VTK mesh.
 * \tparam CELL_NODES The number of nodes in the elements of the VTK mesh.
 */
template <short DIM, short CELL_NODES>
class VtkIO {

private:

    std::vector<std::string> parsed_mesh_;                              /**< The container storing the parsed mesh for random-access to any line */

    std::vector<int> parts_startlines_;                                 /**< The indices of the starting line of the partition sets in the mesh file */

    std::vector<int> nsets_startlines_;                                 /**< The indices of the starting line of the node sets in the mesh file */

    int nodes_startline_;                                               /**< The index of the starting line of the nodes in the parsed mesh file */

    int elems_startline_;                                               /**< The index of the starting line of the elements in the parsed mesh file */

    int elems_type_startline_;                                           /**< The index of the starting line of the elements type in the parsed mesh file */

    CellShape parsed_cells_shape_;                                      /**< The shape of the parsed mesh cells */

    bool parts_exist_;                                                  /**< Conditional to check if the mesh is partitioned */

    bool nsets_exist_;                                                  /**< Conditional to check if the mesh has boundary node sets */


public:

    /**
     * \brief VtkIO constructor.
     */
    inline VtkIO();


    /**
     * \brief VtkIO destructor.
     */
    inline virtual ~VtkIO();


    /**
     * \brief Read a VTK Simple Legacy mesh.
     *
     * The mesh to be readed should be in VTK format (.vtk).
     *
     * \param [in] mesh_filename The filename (full path) of the mesh to be readed.
     * \return [void]
     */
    inline void ReadMeshFrom(const std::string &mesh_filename);


    /**
     * \brief Save mesh in VTK Simple Legacy format.
     * \param [in] mesh The mesh to be saved.
     * \param [in] mesh_filename The filename (full path) where the mesh will be saved.
     * \return [void]
     */
    inline void SaveMesh(const std::vector<Vec<DIM, double> > &nodes, 
                         const std::vector<Cell<DIM, CELL_NODES> > &cells, std::string &mesh_filename);


    /**
     * \brief Load the nodes of the readed mesh in the given nodes' container.
     * \param [out] nodes The container to load the nodes of the mesh.
     * \return [void]
     */
    inline void LoadNodesIn(std::vector<Vec<DIM, double> > &nodes);


    /**
     * \brief Load the elements of the readed mesh in the given elements container.
     * \param [out] tetras The tetrahedra container to load the tetrahedra of the mesh.
     * \return [void]
     */
    inline void LoadElementsIn(std::vector<Cell<DIM, CELL_NODES> > &cells);


    /**
     * \brief Assign the partitions of the readed mesh to the elements of the given container.
     * \param [out] tetras The tetrahedra container where partitions will be assigned to the stored tetrahedra.
     * \return [void]
     */
    //inline void LoadPartitionsIn(std::vector<IMP::Tetrahedron> &tetras);


    /**
     * \brief Load the mesh node sets in the given node sets container.
     * \param [out] node_sets The node sets container where the available mesh node sets will be assigned.
     * \return [void]
     */
    inline void LoadNodeSetsIn(std::unordered_map<std::string, NodeSet> &node_sets);


    /**
     * \brief Conditional to check if the mesh has multiple partitions (domains)
     * \return [bool] True if the mesh has more than one partitions (domains), false differently.
     */
    //inline const bool & PartitionsExist() const { return this->parts_exist_; }


    /**
     * \brief Checks if node sets are available in the mesh.     * 
     * \return [true] Node sets are available in the mesh. 
     * @return [false] There are no node sets in the mesh.
     */
    inline const bool & NodeSetsExist() const { return this->nsets_exist_; }

};


/** \} End of Doxygen Groups*/

} //end of namespace IMP

#endif // IMP_ENGINE_MESH_IO_VTK_IO_HPP_

#include "IMP/engine/mesh_io/vtk_io.tpp"