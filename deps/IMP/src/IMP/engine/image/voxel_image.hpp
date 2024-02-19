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
   \file voxel_image.hpp
   \brief VoxelImage class header file.
   \author Konstantinos A. Mountris
   \date 20/05/2017
*/

#ifndef IMP_ENGINE_IMAGE_VOXEL_IMAGE_HPP_
#define IMP_ENGINE_IMAGE_VOXEL_IMAGE_HPP_

#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/utilities/logger.hpp"
#include "IMP/engine/image/voxel_image_properties.hpp"


#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include <stdexcept>
#include <exception>


namespace IMP {

/** \addtogroup Image \{ */


/**
 * \class VoxelImage
 * \brief Templated class implemmenting a voxelized image.
 */

template <typename VOXELTYPE>
class VoxelImage {

private:
      std::vector<VOXELTYPE> data_;                 /**< The data of the voxelized image */
            
      IMP::Vec<3, int> dimensions_;                 /**< The dimensions of the voxelized image */
      
      IMP::Vec<3, double> spacing_;                 /**< The spacing of the voxelized image */

      IMP::Vec<3, double> origin_;                  /**< The origin of the voxelized image */
      
      int num_voxels_;                              /**< The number of voxels of the voxelized image */
      
      int memory_size_;                             /**< The memory size (number of bytes) of the voxelized image */

      int compressed_bytes_;                        /**< The number of bytes of the voxelized image when is compressed */

      IMP::ImageDataType datatype_;                 /**< The datatype of the voxelized image */
      
      IMP::ImageCompression compression_state_;     /**< The datatype of the voxelized image */


public:     
      
      /**
       * \brief VoxelImage constructor.
       */
      inline VoxelImage();
      
      
      /**
       * \overload
       * \brief VoxelImage copy constructor.
       * \param [in] image The voxelized image to be copied.
       */
      inline VoxelImage(const VoxelImage &image);


      /**
       * \overload
       * \brief VoxelImage move constructor.
       * \param [in] image The voxelized image to be moved.
       */
      inline VoxelImage(VoxelImage &&image);

      
      /**
       * \brief VoxelImage destructor.
       */
      inline virtual ~VoxelImage();

      
      /**
       * \brief Initialize an empty voxelized image.
       * \param [in] dimensions The dimensions of the voxelized image.
       * \param [in] spacing The voxel spacing of the voxelized image.
       * \param [in] origin The origin of the voxelized image.
       * \return [void]
       */
      inline void InitEmptyImage(const IMP::Vec<3, int> &dimensions,
                          const IMP::Vec<3, double> &spacing,
                          const IMP::Vec<3, double> &origin);
      
      
      /**
       * \brief Set dimensions of the voxelized image.
       * \param [in] x_dim The X dimension of the voxelized image.
       * \param [in] y_dim The Y dimension of the voxelized image.
       * \param [in] z_dim The Z dimension of the voxelized image.
       * \return [void]
       * \throws std::invalid_argument
       */
      inline void SetDimensions(const int &x_dim, const int &y_dim, const int &z_dim);
      
      
      /**
       * \brief Set spacing of the voxelized image.
       * \param [in] x_spacing The X spacing of the voxelized image.
       * \param [in] y_spacing The Y spacing of the voxelized image.
       * \param [in] z_spacing The Z spacing of the voxelized image.
       * \return [void]
       * \throws std::invalid_argument
       */
      inline void SetSpacing(const double &x_spacing, const double &y_spacing, const double &z_spacing);


      /**
       * \brief Set origin of the voxelized image.
       * \param [in] x_origin The X origin of the voxelized image.
       * \param [in] y_origin The Y origin of the voxelized image.
       * \param [in] z_origin The Z origin of the voxelized image.
       * \return [void]
       */
      inline void SetOrigin(const double &x_origin, const double &y_origin, const double &z_origin);


      /**
       * \brief Set the number of voxels of the voxelized image.
       * \param [in] num_voxels The number of voxels of the voxelized image.
       * \return [void]
       * \throws std::invalid_argument
       */
      inline void SetNumVoxels(const int &num_voxels);


      /**
       * \brief Set the memory size occupied by the voxelized image.
       * \param [in] memory_size The memory size occupied by the voxelized image.
       * \return [void]
       * \throws std::invalid_argument
       */
      inline void SetMemorySize(const int &memory_size);


      /**
       * \brief Set the datatype of the voxelized image.
       * \param [in] datatype The datatype of the voxelized image.
       * \return [void]
       */
      inline void SetDataType(const IMP::ImageDataType &datatype);


      /**
       * \brief Set the compression state of the voxelized image.
       * \param [in] compress_state The compression state of the voxelized image.
       * \return [void]
       */
      inline void SetCompression(const IMP::ImageCompression &compress_state);


      /**
       * \brief Set the number of bytes of the voxelized image when is compressed.
       * \param [in] compressed_bytes The number of bytes of the voxelized image when is compressed.
       * \return [void]
       * \throws std::invalid_argument
       */
      inline void SetCompressedBytes(const int &compressed_bytes);


      /**
       * \brief Check if the image is initialized.
       *
       * The image is considered initialized if the number of voxels is gretear than zero,
       * it is equal to the product of the three image dimensions and the memory size of the image
       * is equal to: [number of voxels x sizeof(VOXELTYPE)], where VOXELTYPE is the data type of the image.
       * \return [bool] True if the image is initialized.
       */
      inline bool IsInitialized() const;
      
      
      /**
       * \brief Load voxelized image from file.
       *
       * Supported formats: [.h33 | .hdr | .inr | .mha | .mhd | .nii].
       *
       * \param[in] img_filename The filename (with full path) of the image to load.
       */
      inline void LoadFrom(const std::string &filename);


      /**
       * \brief Save voxelized image to file.
       *
       * Supported formats: [.h33 | .hdr | .inr | .mha | .mhd | .nii].
       *
       * \param[in] img_filename The filename (with full path) where the image will be saved.
       * \param[in] compression The desired compression state (ON/OFF) for the output image. [Default: OFF].
       */
      inline void SaveTo(const std::string &img_filename, ImageCompression compression = ImageCompression::OFF);


      /**
       * \brief Write access to the data of the voxelized image.
       * \return [std::vector<VOXELTYPE>] The data of the voxelized image with write access.
       */
      inline std::vector<VOXELTYPE> & EditData() { return this->data_; }
      
      
      /**
       * \brief Read-only access to the data of the voxelized image.
       * \return [std::vector<VOXELTYPE>] The data of the voxelized image with read-only access.
       */
      inline const std::vector<VOXELTYPE> & Data() const { return this->data_; }
      
      
      /**
       * \brief Get the voxelized image dimensions.
       * \return [IMP::Vec<3, int>] The the voxelized image dimensions.
       */
      inline const IMP::Vec<3, int> & Dimensions() const { return this->dimensions_; }
      
      
      /**
       * \brief Get the voxel spacing of the voxelized image in the three dimensions.
       * \return [IMP::Vec<3, double>] The voxel spacing of the voxelized image in the three dimensions.
       */
      inline const IMP::Vec<3, double> & Spacing() const { return this->spacing_; }
      
      
      /**
       * \brief Get the origin from the origin of the voxelized image in the three dimensions.
       * \return [IMP::Vec<3, double>] The origin from the origin of the voxelized image in the three dimensions.
       */
      inline const IMP::Vec<3, double> & Origin() const { return this->origin_; }
      

      /**
       * \brief Get the voxel volume of the voxelized image.
       * \return [double] The voxel volume of the voxelized image.
       */
      inline const double & VoxelVolume() const { return this->spacing_[0]*this->spacing_[1]*this->spacing_[2]; }
   

      /**
       * \brief Get the number of voxels of the voxelized image.
       * \return [int] The number of the voxels of the voxelized image.
       */
      inline const int  & NumVoxels() const { return this->num_voxels_; }


      /**
       * \brief Get the size in memory (number of bytes) of the voxelized image.
       * \return [int] The memory size (number of bytes) of the voxelized image.
       */
      inline const int & MemorySize() const { return this->memory_size_; }


      /**
       * \brief Get the datatype of the voxelized image.
       * \return [IMP::ImageDataType] The datatype of the voxelized image.
       */
      inline const IMP::ImageDataType & DataType() const { return this->datatype_; }


      /**
       * \brief Get the compression state of the voxelized image.
       * \return [IMP::ImageCompression] The compression state of the voxelized image.
       */
      inline const IMP::ImageCompression & Compression() const { return this->compression_state_; }


      /**
       * \brief Get the number of bytes of the voxelized image when is compressed.
       * \return [int] The number of bytes of the voxelized image when is compressed.
       */
      inline const int & CompressedBytes() const { return this->compressed_bytes_; }


      /**
       * \brief Get a pointer to the first voxel of the voxelized image.
       * \return [VOXELTYPE pointer] The pointer to the first voxel of the voxelized image.
       */
      inline VOXELTYPE *DataPtr() { return &(*this->data_.begin()); }


      /**
       * \brief Equal to operator.
       *
       * Compares voxelized images for equality.
       *
       * \param [in] image The voxelized image to compare.
       * \return [bool] TRUE if voxelized images are identical.
       */
      inline bool operator == (const VoxelImage &image) const;


      /**
       * \brief Not equal to operator.
       *
       * Compares voxelized images for inequality.
       *
       * \param [in] image The voxelized image to compare.
       * \return [bool] TRUE if voxelized images are not identical.
       */
      inline bool operator != (const VoxelImage &image) const;


      /**
       * \brief Copy assignment operator.
       *
       * Assigns all the properties of a given voxelized image (data, dimensions, spacing, origin, etc.) by copying.
       *
       * \param [in] image The voxelize image to assign.
       * \return [IMP::VoxelImage] The assigned voxelized image.
       */
      inline VoxelImage<VOXELTYPE> & operator = (const VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Move assignment operator.
       *
       * Assigns all the properties of a given voxelized image (data, dimensions, spacing, origin, etc.) by moving.
       *
       * \param [in] image The voxelize image to assign.
       * \return [IMP::VoxelImage] The assigned voxelized image.
       */
      inline VoxelImage<VOXELTYPE> & operator = (VoxelImage<VOXELTYPE> &&image);


      inline VOXELTYPE operator [] (const int &id) const;


      inline VOXELTYPE operator () (const int &i, const int &j, const int &k) const;
      
      
      /**
       * \brief Clear the data and characheristics of the voxelized image.
       * \return [void]
       */
      inline void Clear();

};


/** \} End of Doxygen Groups*/
} //end of namespace IMP

#endif // IMP_ENGINE_IMAGE_VOXEL_IMAGE_HPP_

#include "IMP/engine/image/voxel_image.tpp"