/*
 * CLOUDEA - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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



/*!
   \file febio_io.hpp
   \brief Febio input/output class header file.
   \author Konstantinos A. Mountris
   \date 08/01/2018
*/

#ifndef CLOUDEA_MESH_IO_FEBIO_IO_HPP_
#define CLOUDEA_MESH_IO_FEBIO_IO_HPP_


#include "CLOUDEA/engine/vectors/vec3.hpp"
#include "CLOUDEA/engine/elements/element_properties.hpp"
#include "CLOUDEA/engine/elements/node.hpp"
#include "CLOUDEA/engine/elements/tetrahedron.hpp"
#include "CLOUDEA/engine/sets/node_set.hpp"
#include "CLOUDEA/engine/mesh/mesh_properties.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"
#include "CLOUDEA/engine/utilities/tinyxml2.h"

// #include <tinyxml2.h>

#include <vector>
#include <map>

#include <utility>
#include <cctype>

#include <exception>
#include <stdexcept>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <algorithm>
#include <string>

namespace CLOUDEA {

/*!
 *  \addtogroup MeshIO
 *  @{
 */

/*!
 * \class FebioIO
 * \brief Class implemmenting input/output functionality for mesh in FeBio format (.feb).
 *
 */
class FebioIO {
public:

    /*!
     * \brief FebioIO default constructor.
     */
    FebioIO();


    /*!
     * \brief FebioIO default destructor.
     */
    virtual ~FebioIO();


    /*!
     * \brief Load a FeBio mesh.
     *
     * The mesh to be loaded should be in FeBio format (.feb).
     *
     * \param [in] mesh_filename The filename (full path) of the mesh to be loaded.
     * \return [void]
     */
    void LoadMeshFrom(const std::string &mesh_filename);


    /*!
     * \brief Load the vertices of the loaded mesh in the given vertices container.
     * \param [out] vertices The vertices container to load the vertices of the mesh.
     * \return [void] The vertices of the loaded mesh.
     */
    void LoadNodesIn(std::vector<Node> &nodes);


    /*!
     * \brief Load the elements of the loaded mesh in the given elements container.
     * \param [out] tetras The tetrahedra container to load the tetrahedra of the mesh.
     * \return [void]
     */
    void LoadElementsIn(std::vector<Tetrahedron> &tetras);


    /*!
     * \brief Assign the boundary node sets of the loaded mesh to the given container.
     * \param node_sets The container where the loaded mesh nodesets will be stored.
     * \return [void]
     */
    void LoadBoundarySetsIn(std::vector<NodeSet> &node_sets);


    /*!
     * \brief Conditional to check if the mesh has multiple boundaries.
     * \return [bool] True if the mesh has more than one boundaries, false differently.
     */
    inline const bool & BoundariesExist() const { return this->boundaries_exist; }



private:

    tinyxml2::XMLDocument input_mesh_;                       /*!< The parsed mesh. */

    std::vector<std::pair<int,int> > offsetted_nodes_;       /*!< The indices of the offsetted nodes from the storage order and their offsets. */

    std::vector<std::pair<int,int> > offsetted_elems_;       /*!< The indices of the offsetted elements from the storage order and their offsets. */

    bool is_mesh_loaded;

    bool boundaries_exist;                                   /*!< The conditional of the mesh boundaries' existence */

};


/*! @} End of Doxygen Groups*/

} //end of namespace CLOUDEA


#endif // CLOUDEA_MESH_IO_FEBIO_IO_HPP_
