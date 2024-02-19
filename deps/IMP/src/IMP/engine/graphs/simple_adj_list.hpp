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
   \file simple_adj_list.hpp
   \brief SimpleAdjList class header file.
   \author Konstantinos A. Mountris
   \date 24/11/2020
*/

#ifndef IMP_ENGINE_GRAPHS_SIMPLE_ADJ_LIST_HPP_
#define IMP_ENGINE_GRAPHS_SIMPLE_ADJ_LIST_HPP_

#include <vector>

namespace IMP {

/** \addtogroup Graphs \{ */


/**
 * \class SimpleAdjList
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting an undirected adjacency list graph.
 */
class SimpleAdjList {

private:

    std::vector<std::vector<int>> links_;       /**< The undirected adjacency list of linked objects */

    int links_num_;                             /**< The number of listed links in the ajacency list */


public:

    /**
     * \brief The SimpleAdjList constructor.
     */
    inline SimpleAdjList() : links_(), links_num_(0) {}


    /**
     * \brief The SimpleAdjList constructor.
     */
    inline SimpleAdjList(int links_num) : links_(), links_num_(links_num) { links_.resize(links_num); }

    
    /**
     * \brief The SimpleAdjList constructor.
     */
    inline SimpleAdjList(const SimpleAdjList &adj_list) : links_(adj_list.links_), links_num_(adj_list.links_num_) {}


    /**
     * \brief The SimpleAdjList destructor.
     */
    inline virtual ~SimpleAdjList() {}


    /**
     * @brief Set the number of listed links.
     * 
     * @param num 
     */
    void SetLinksNum(std::size_t num);


    /**
     * @brief Add a connection in a listed link.
     * 
     * @param i The index of the link.
     * @param j The index of the linked object in the link.
     */
    void AddConnection(std::size_t i, std::size_t j);


    SimpleAdjList & operator = (const SimpleAdjList &adj_list);


    /**
     * @brief Get the listed links. 
     * @return const std::vector<std::vector<int>>& 
     */
    inline const std::vector<std::vector<int>> & Links() const { return this->links_; }


    /**
     * @brief Get the listed link with the specified index. 
     * Fast access without range check.
     * @return const std::vector<std::vector<int>>& 
     */
    inline const std::vector<int> & Links(std::size_t id) const { return this->links_[id]; }


    /**
     * @brief Get the listed link with the specified index. 
     * Slow access with range check.
     * @return const std::vector<std::vector<int>>& 
     */
    inline const std::vector<int> & LinksAt(std::size_t id) const { return this->links_.at(id); }


    /**
     * @brief Get the number of the linked lists.
     * 
     * @return int 
     */
    inline int LinksNum() const { return this->links_num_; }

};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_GRAPHS_SIMPLE_ADJ_LIST_HPP_ 
