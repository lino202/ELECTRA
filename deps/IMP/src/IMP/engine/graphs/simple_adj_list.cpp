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


#include "IMP/engine/graphs/simple_adj_list.hpp"

namespace IMP {


void SimpleAdjList::SetLinksNum(std::size_t num)
{
    this->links_num_ = num;
    this->links_.clear();
    this->links_.resize(num);
    this->links_.shrink_to_fit();
}


void SimpleAdjList::AddConnection(std::size_t i, std::size_t j)
{
    this->links_[i].emplace_back(j);
}


SimpleAdjList & SimpleAdjList::operator = (const SimpleAdjList &adj_list)
{
    this->links_ = adj_list.links_;
    this->links_num_ = adj_list.links_num_;
    return *this;
}


} // End of namespace IMP.
