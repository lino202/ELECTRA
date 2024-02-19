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
   \file mhd_io.hpp
   \brief MhdIo class header file.
   \author Konstantinos A. Mountris
   \date 23/05/2017
*/

#ifndef IMP_ENGINE_IMAGE_IO_MHD_IO_HPP_
#define IMP_ENGINE_IMAGE_IO_MHD_IO_HPP_

#include "IMP/engine/vectors/vec.hpp"
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
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/copy.hpp>

namespace IMP {

/** \addtogroup ImageIO \{ */


/**
 * \class MhdIO
 * \brief Class implemmenting input/output for voxelized images in MetaImage format [.mhd].
 */


class MhdIO{
public:

      /**
       * \brief MhdIO constructor.
       */
      inline MhdIO();


      /**
       * \brief MhdIO destructor.
       */
      inline virtual ~MhdIO();


      /**
       * \brief Load voxelized image.
       *
       * The voxelized image format must be .mhd.
       *
       * \param[in] filename The filename (with full path) of the image to load.
       * \param[out] image The voxelized image where image data will be loaded.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void Load(const std::string &filename, VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Save voxelized image to MetaImage file format (.mhd).
       * \param[in] image The voxelized image to be saved.
       * \param[in] filename The filename (with full path) to save the image.
       * \param[in] compression The desired compression state (ON/OFF) for the output image.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression);


      /**
       * \brief Read binary (.raw) data of the voxelized image.
       * \param[in] raw_filename The filename (with full path) of the binary (.raw) file.
       * \param[out] image The voxelized image where raw data will be stored.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void ReadRawData(const std::string &raw_filename, VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Read the compressed binary (.zraw) data of the voxelized image.
       * \param[in] zraw_filename The filename (with full path) of the compressed binary (.zraw) file.
       * \param[out] image The voxelized image where raw data will be stored after decompression.
       * \return [void]
       */
      template <typename VOXELTYPE>
      inline void ReadZrawData(const std::string &zraw_filename, VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Write binary (.raw) data of the voxelized image.
       * \param[in] image The voxelized image to write in the binary output file.
       * \param[in] raw_filename The filename (with full path) of the birary (.raw) output file.
       * \return [void]
       */
      template<typename VOXELTYPE>
      inline void WriteRawData(VoxelImage<VOXELTYPE> &image, const std::string &raw_filename);


      /**
       * \brief Write compressed binary (.zraw) data of the voxelized image.
       * \param[in] image The voxelized image to write in the compressed birary output file.
       * \param[in] zraw_filename The filename (with full path) of the compressed birary (.zraw) output file.
       * \return [void]
       */
      template<typename VOXELTYPE>
      inline void WriteZRawData(VoxelImage<VOXELTYPE> &image, const std::string &zraw_filename);



private:

      /**
       * \brief Extract a string from a MetaImage header (.mhd).
       * \param[in] line The line of the MetaImage header (.mhd) to be processed.
       * \param[out] attribute The part of the string to be extracted.
       * \return [void]
       */
      inline void ExtractStringFromHeader(const std::string &line, std::string &attribute);


      /**
       * \brief Extract a 3-dimensional attribute from a MetaImage header (.mhd).
       * \param[in] line The line of the MetaImage header (.mhd) to be processed.
       * \param[out] vec The IMP::Vec3 that will be extracted.
       * \return [void]
       */
      template<typename DATATYPE>
      inline void ExtractVecFromHeader(const std::string &line, Vec<3, DATATYPE> &vec);


};


/** @} End of Doxygen Groups*/
} //end of namespace IMP

#endif //IMP_ENGINE_IMAGE_IO_MHD_IO_HPP_

#include "IMP/engine/image_io/mhd_io.tpp"
