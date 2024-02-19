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
   \file cell_set.hpp
   \brief CellSet class header file.
   \author Konstantinos A. Mountris
   \date 18/09/2019
*/


#ifndef IMP_ENGINE_TOPOLOGY_CELL_SET_HPP_
#define IMP_ENGINE_TOPOLOGY_CELL_SET_HPP_

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

namespace IMP {

/** \addtogroup Topology-Sets \{ */


/**
 * \class CellSet
 * \brief Class implemmenting a set of cell indices and related attributes to the set.
 */
class CellSet {

private:

    std::vector<int> cell_ids_;      /**< The indices of the cells belonging to the set */

    std::string name_;               /**< The name attribute of the cell set */

public:

    /**
     * \brief The CellSet default constructor.
     */
    CellSet();


    /**
     * \brief The CellSet copy constructor.
     * \param [in] cell_set The node set to be copied.
     */
    CellSet(const CellSet &cell_set);


    /**
     * \brief The CellSet destructor.
     */
    virtual ~CellSet();


    /**
     * \brief Set the cell set's name attribute.
     * \param [in] name The name attribute of the cell set. 
     * \return [void]
     */
    inline void SetName(const std::string &name) { this->name_ = name; }


    /**
     * \brief Clear the cell set. 
     * \return [void]
     */
    inline void Clear() { this->name_ = ""; this->cell_ids_.clear(); }


    /**
     * \brief Get the cell indices of the cell set with editing access. 
     * \return [std::vector<int>&] The cell set indices with editing access. 
     */
    inline std::vector<int> & EditCellIds() { return this->cell_ids_; }


    /**
     * \brief Get the cell indices of the cell set.
     * \return [const std::vector<int>&] The cell indices of the cell set.
     */
    inline const std::vector<int> & CellIds() const { return this->cell_ids_; }


    /**
     * \brief Get the name attribute of the cell set.
     * @return [const std::string&] The name attribute of the cell set. 
     */
    inline const std::string & Name() const { return this->name_; }


    /**
     * \brief Equal to operator.
     *
     * Compares cell sets for equality.
     *
     * \param [in] cell_set The cell set to compare.
     * \return [bool] TRUE: cell sets are identical | FALSE: otherwise.
     */
    bool operator == (const CellSet &cell_set) const;


    /**
     * \brief Not equal to operator.
     *
     * Compares cell sets for inequality.
     *
     * \param [in] cell_set The cell set to compare.
     * \return [bool] TRUE: cell sets are "not" identical | FALSE: otherwise.
     */
    bool operator != (const CellSet &cell_set) const;


    /**
     * \brief Assignment operator.
     *
     * Assigns all the properties of a given cell set (cell set's name, cell indices).
     *
     * \param [in] cell_set The cell set to assign.
     * \return [IMP::CellSet&] The assigned cell set.
     */
    CellSet & operator = (const CellSet &cell_set);

};


/** \} End of Doxygen Groups*/
} //namespace IMP

#endif //IMP_ENGINE_TOPOLOGY_CELL_SET_HPP_