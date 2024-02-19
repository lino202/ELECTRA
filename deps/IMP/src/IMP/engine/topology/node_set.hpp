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
   \file node_set.hpp
   \brief NodeSet class header file.
   \author Konstantinos A. Mountris
   \date 21/03/2019
*/


#ifndef IMP_ENGINE_TOPOLOGY_NODE_SET_HPP_
#define IMP_ENGINE_TOPOLOGY_NODE_SET_HPP_

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

namespace IMP {

/** \addtogroup Topology-Sets \{ */


/**
 * \class NodeSet
 * \brief Class implemmenting a set of node indices and related attributes to the set.
 */
class NodeSet {

private:

    std::vector<int> node_ids_;      /**< The indices of the nodes belonging to the set */

    std::string name_;               /**< The name attribute of the node set */


public:

    /**
     * \brief The NodeSet default constructor.
     */
    NodeSet();


    /**
     * \brief The NodeSet copy constructor.
     * \param [in] node_set The node set to be copied.
     */
    NodeSet(const NodeSet &node_set);


    /**
     * \brief The NodeSet destructor.
     */
    virtual ~NodeSet();


    /**
     * \brief Set the node set's name attribute.
     * \param [in] name The name attribute of the node set.
     * \return [void]
     */
    inline void SetName(const std::string &name) { this->name_ = name; }


    /**
     * \brief Set the node set's name attribute and node indices.
     * \param [in] name The name of the node set.
     * \param [in] node_ids The node indices of the nodeset.
     */
    inline void Set(const std::string &name, const std::vector<int> &node_ids) { this->name_ = name; this->node_ids_ = node_ids; }


    /**
     * \brief Clear the node set.
     * \return [void]
     */
    inline void Clear() { this->name_ = ""; this->node_ids_.clear(); }


    /**
     * \brief Get the node indices of the node set with editing access.
     * \return [std::vector<int>&] The node set indices with editing access.
     */
    inline std::vector<int> & EditNodeIds() { return this->node_ids_; }


    /**
     * \brief Get the node indices of the node set.
     * \return [const std::vector<int>&] The node indices of the node set.
     */
    inline const std::vector<int> & NodeIds() const { return this->node_ids_; }


    /**
     * \brief Get the name attribute of the node set.
     * @return [const std::string&] The name attribute of the node set.
     */
    inline const std::string & Name() const { return this->name_; }


    /**
     * \brief Equal to operator.
     *
     * Compares node sets for equality.
     *
     * \param [in] node_set The node set to compare.
     * \return [bool] TRUE: node sets are identical | FALSE: otherwise.
     */
    bool operator == (const NodeSet &node_set) const;


    /**
     * \brief Not equal to operator.
     *
     * Compares node sets for inequality.
     *
     * \param [in] node_set The node set to compare.
     * \return [bool] TRUE: node sets are "not" identical | FALSE: otherwise.
     */
    bool operator != (const NodeSet &node_set) const;


    /**
     * \brief Assignment operator.
     *
     * Assigns all the properties of a given node set (node set's name, node indices).
     *
     * \param [in] node_set The node set to assign.
     * \return [IMP::NodeSet&] The assigned node set.
     */
    NodeSet & operator = (const NodeSet &node_set);

};


/** \} End of Doxygen Groups*/
} //namespace IMP

#endif //IMP_ENGINE_TOPOLOGY_NODE_SET_HPP_
