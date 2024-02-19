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
   \file cell_props.hpp
   \brief Collection of cell properties header file.
   \author Konstantinos A. Mountris
   \date 30/05/2019
*/

#ifndef IMP_ENGINE_ELEMENTS_CELL_PROPS_HPP_
#define IMP_ENGINE_ELEMENTS_CELL_PROPS_HPP_


namespace IMP {

/** \addtogroup Elements \{ */

/**
 * \enum
 */
enum struct CellShape: int { unknown = 0,    /**< Cell type -> unknown */
                             mixed = 1,      /**< Cell type -> mixed */
                             node = 2,       /**< Cell type -> node */
                             edge = 3,       /**< Cell type -> edge */
                             tri = 4,        /**< Cell type -> triangle */
                             quad = 5,       /**< Cell type -> quadrilateral */
                             tet = 6,        /**< Cell type -> tetrahedron */
                             hex = 7,        /**< Cell type -> hexahedron */
                             poly = 8        /**< Cell type -> polygonal/polyhedral */       
                           };
/** \} End of Doxygen Groups*/

} //end of namespace IMP

#endif //IMP_ELEMENTS_CELL_PROPS_HPP_
