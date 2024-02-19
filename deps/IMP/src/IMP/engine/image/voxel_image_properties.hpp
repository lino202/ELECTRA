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
   \file voxel_image_properties.hpp
   \brief Collection of voxelized image properties header file.
   \author Konstantinos A. Mountris
   \date 12/06/2017
*/

#ifndef IMP_ENGINE_IMAGE_VOXEL_IMAGE_PROPERTIES_HPP_
#define IMP_ENGINE_IMAGE_VOXEL_IMAGE_PROPERTIES_HPP_


namespace IMP {

/** \addtogroup Image \{ */


/**
 * \enum ImageDataType
 * \brief Enumeration for various datatypes that are supported by the VoxelImage class.
 */
enum struct ImageDataType: int {no_type = 0,                    /**< Image data type -> no type. */
                                char_type = 1,                  /**< Image data type -> signed char. */
                                uchar_type = 2,                 /**< Image data type -> unsigned char. */
                                short_int_type = 3,             /**< Image data type -> signed short int. */
                                ushort_int_type = 4,            /**< Image data type -> unsigned short int. */
                                int_type = 5,                   /**< Image data type -> signed int. */
                                uint_type = 6,                  /**< Image data type -> unsigned int. */
                                long_int_type = 7,              /**< Image data type -> long int. */
                                ulong_int_type = 8,             /**< Image data type -> unsigned long int. */
                                float_type = 9,                 /**< Image data type -> float. */
                                double_type = 10,               /**< Image data type -> double. */
                                bool_type = 11                  /**< Image data type -> bool. */
                               };


/**
 * \enum ImageCompressionType
 * \brief Enumeration for the compression state of the voxelized image implemented by VoxelImage.
 */
enum struct ImageCompression: bool {OFF = false,   /**< Image Compression type -> OFF when no compression is required. */
                                    ON = true      /**< Image Compression type -> ON when compression is required. */
                                   };



/** \} End of Doxygen Groups*/

} //end of namespace IMP

#endif //IMP_ENGINE_IMAGE_VOXEL_IMAGE_PROPERTIES_HPP_
