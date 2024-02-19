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
   \file h33_io.hpp
   \brief H33Io class header file.
   \author Konstantinos A. Mountris
   \date 28/06/2017
*/

#ifndef IMP_ENGINE_IMAGE_IO_H33_IO_HPP_
#define IMP_ENGINE_IMAGE_IO_H33_IO_HPP_

#include "IMP/engine/image/voxel_image.hpp"

#include <string>
#include <stdexcept>
#include <exception>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstddef>

#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/counter.hpp>

namespace IMP {

/** \addtogroup ImageIO \{ */


/**
 * \class H33IO
 * \brief Class implemmenting input/output for voxelized images in interfile format [.h33].
 */

/**
 * \todo Compression read write implementation
*/
class H33IO{
public:

      /**
       * \brief H33IO constructor.
       */
      inline H33IO();


      /**
       * \brief H33IO destructor.
       */
      inline virtual ~H33IO();


      /**
       * \brief Load voxelized image.
       *
       * The voxelized image format must be .h33 (interfile).
       *
       * \param[in] filename The filename (with full path) of the image to load.
       * \param[out] image The voxelized image where image data will be loaded.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void Load(const std::string &filename, VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Save voxelized image to Interfile file format (.h33).
       * \param[in] image The voxelized image to be saved.
       * \param[in] filename The filename (with full path) to save the image.
       * \param[in] compression The desired compression state (ON/OFF) for the output image.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression);


      /**
       * \brief Read binary data (.i33) of the voxelized image.
       * \param[in] i33_filename The filename (with full path) of the binary (.i33) file.
       * \param[out] image The voxelized image where .i33 data will be stored.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void ReadI33Data(const std::string &i33_filename, VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Read compressed binary data (.i33.gz) of the voxelized image.
       *
       * Supported compressing by the gzip library.
       *
       * \param[in] i33_gz_filename The filename (with full path) of the compressed binary data (.i33.gz).
       * \param[out] image The voxelized image where compressed binary data (.i33.gz) will be stored.
       * \return [void]
       */
      // template <typename VOXELTYPE>
      // inline void ReadI33GzData(const std::string &i33_gz_filename, VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Write binary (.i33) data of the voxelized image.
       * \param[in] image The voxelized image to write in output file.
       * \param[in] i33_filename The filename (with full path) of the birary (.i33) output file.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void WriteI33Data(VoxelImage<VOXELTYPE> &image, const std::string &i33_filename);


      /**
       * \brief Write compressed binary (.i33.gz) data of the voxelized image.
       * \param[in] image The voxelized image to write in the compressed birary output file.
       * \param[in] i33_gz_filename The filename (with full path) of the compressed birary (.i33.gz) output file.
       * \return [void]
       */
      // template <typename VOXELTYPE>
      // inline void WriteI33GzData(VoxelImage<VOXELTYPE> &image, const std::string &i33_gz_filename);




private:

      /**
       * \brief Extract an attribute from a Inria image header (.h33).
       * \param[in] line The line of the Interfile image header (.h33) to be processed.
       * \param[out] attribute The attribute that will be extracted.
       * \return [void]
       * \return [void]
       */
      template<typename DATATYPE>
      inline void ExtractAttributeFromHeader(const std::string &line, DATATYPE &attribute);


};


/** \} End of Doxygen Groups*/
} //end of namespace IMP

#endif //IMP_ENGINE_IMAGE_IO_H33_IO_HPP_

#include "IMP/engine/image_io/h33_io.tpp"
