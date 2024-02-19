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

#include "IMP/engine/topology/cell_set.hpp"


namespace IMP {


CellSet::CellSet() : cell_ids_(), name_("")
{}


CellSet::CellSet(const CellSet &cell_set) : cell_ids_(), name_("")
{
    this->name_ = cell_set.name_;
    this->cell_ids_ = cell_set.cell_ids_;
}


CellSet::~CellSet()
{}


bool CellSet::operator == (const CellSet &cell_set) const
{
    // Compare cell sets for equality.
    return ((this->name_ == cell_set.name_) &&
            (this->cell_ids_ == cell_set.cell_ids_)
           );
}


bool CellSet::operator != (const CellSet &cell_set) const
{
    // Compare cell sets for inequality.
    return !(*this == cell_set);
}


CellSet & CellSet::operator = (const CellSet &cell_set)
{
    if (this != &cell_set) {
        // Assign values from cell_set.
        this->name_ = cell_set.name_;
        this->cell_ids_ = cell_set.cell_ids_;
    }

    return *this;
}


} // End of namespace IMP