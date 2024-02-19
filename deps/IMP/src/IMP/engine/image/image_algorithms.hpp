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
   \file image_algorithms.hpp
   \brief Collection of algorithms acting on images.
   \author Konstantinos A. Mountris
   \date 25/06/2018
*/

#ifndef IMP_ENGINE_IMAGE_IMAGE_ALGORITHMS_HPP_
#define IMP_ENGINE_IMAGE_IMAGE_ALGORITHMS_HPP_

#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/utilities/logger.hpp"
#include "IMP/engine/image/voxel_image_properties.hpp"
#include "IMP/engine/image/voxel_image.hpp"

// #include <CGAL/ImageIO.h>
// #include <CGAL/Image_3.h>

#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <stdexcept>
#include <exception>


namespace IMP {


/** \addtogroup Image \{ */


/**
 * \brief Convert a VoxelImage to CGAL::Image_3.
 * \param [in] vox_image The IMP::VoxelImage to convert to CGAL::Image_3.
 * \return [CGAL::Image_3] The converted CGAL::Image_3.
 */
// template <typename VOXELTYPE>
// inline CGAL::Image_3 CgalImage_3FromIMPVoxelImage(const VoxelImage<VOXELTYPE> &vox_image);


/** \} End of Doxygen Groups*/
} // End of namespace IMP

#endif //IMP_ENGINE_IMAGE_IMAGE_ALGORITHMS_HPP_

#include "IMP/engine/image/image_algorithms.tpp"
