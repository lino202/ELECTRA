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
   \file nii_io.hpp
   \brief NiiIO class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2017
*/

#ifndef IMP_ENGINE_IMAGE_IO_NII_IO_HPP_
#define IMP_ENGINE_IMAGE_IO_NII_IO_HPP_

#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/image/voxel_image.hpp"
#include "IMP/engine/utilities/logger.hpp"
#include "IMP/engine/utilities/timer.hpp"

#include <string>
#include <stdexcept>
#include <exception>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <cstdint>

#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/counter.hpp>
#include <boost/iostreams/copy.hpp>

namespace IMP {

/** \addtogroup ImageIO \{ */


/**
 * \class NiiIO
 * \brief Class implemmenting input/output for voxelized images in NifTI-1 and NifTI-2 formats [.hdr | .nii].
 */


class NiiIO{
public:

      /**
       * \brief NiiIO constructor.
       */
      inline NiiIO();


      /**
       * \brief NiiIO destructor.
       */
      inline virtual ~NiiIO();


      /**
       * \brief Load voxelized image.
       *
       * The voxelized image format must be in NifTI-1 or NifTI-2 format [.hdr | .nii].
       *
       * \param[in] filename The filename (with full path) of the image to load.
       * \param[out] image The voxelized image where image data will be loaded.
       * \return [void]
       *
       * \todo Read gzip compressed files.
       * \bug Is reading correct gzipped img file??
       */
      template <typename VOXELTYPE>
      inline void Load(const std::string &filename, VoxelImage<VOXELTYPE> &image);


      /**
       * \brief Save voxelized image to Inria file format (.mhd).
       * \param[in] image The voxelized image to be saved.
       * \param[in] filename The filename (with full path) to save the image.
       * \param[in] compression The desired compression state (ON/OFF) for the output image.
       * \return [void]
       *
       * \note Currently output in NifTI-1 format is implemented.
       */
      template <typename VOXELTYPE>
      inline void Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression);

};


/**
 * \struct Nifti1Header
 * \brief Structure collecting the attributes of the 348 bytes header of a voxelized image in NifTI-1 format [.hdr | .nii].
 */

typedef struct Nifti1Header {

    /**
     * \brief The Nifti1Header constructor.
     */
    inline Nifti1Header() : hdr_size_(0), data_type_{'\0'}, db_name_{'\0'}, extents_(0), session_error_(0), regular_('\0'),
        dim_info_('\0'), dim_{0}, intent_p1_(0.f), intent_p2_(0.f), intent_p3_(0.f), intent_code_(0), datatype_(0),
        bitpix_(0), slice_start_(0), pixdim_{0.f}, vox_offset_(0.f), scl_slope_(0.f), scl_inter_(0.f), slice_end_(0),
        slice_code_('\0'), xyzt_units_('\0'), cal_max_(0.f), cal_min_(0.f), slice_duration_(0.f), toffset_(0.f), glmax_(0),
        glmin_(0), descrip_{'\0'}, aux_file_{'\0'}, qform_code_(0), sform_code_(0), quatern_b_(0.f), quatern_c_(0.f),
        quatern_d_(0.f), qoffset_x_(0.f), qoffset_y_(0.f), qoffset_z_(0.f), srow_x_{0.f}, srow_y_{0.f}, srow_z_{0.},
        intent_name_{'\0'}, magic_{'\0'}
    {}


    int   hdr_size_;            /**< Size of the header in bytes. It must be 348. */

    char  data_type_[10];       /**< Unused variable. */

    char  db_name_[18];         /**< Unused variable. */

    int   extents_;             /**< Unused variable. */

    short session_error_;       /**< Unused variable. */

    char  regular_;             /**< Unused variable. */

    char  dim_info_;            /**< The ordering of MRI slices. */

    short dim_[8];              /**< The image dimensions.*/

    float intent_p1_;           /**< The 1st intent parameter. */

    float intent_p2_;           /**< The 2nd intent parameter. */

    float intent_p3_;           /**< The 3rd intent parameter. */

    short intent_code_;         /**< The NIFTI_INTENT_* code. */

    short datatype_;            /**< The type of the stored data. */

    short bitpix_;              /**< The number of bits per voxel. */

    short slice_start_;         /**< The index of the first slice. */

    float pixdim_[8];           /**< The image spacings. */

    float vox_offset_;          /**< The offset into the .nii file where the image data start. */

    float scl_slope_;           /**< The slope of the data scaling. */

    float scl_inter_;           /**< The offset of the data scaling. */

    short slice_end_;           /**< The index of the last slice. */

    char  slice_code_;          /**< The slice timing order. */

    char  xyzt_units_;          /**< The units of the image spacings -> pixdim[1..4]. */

    float cal_max_;             /**< The maximum display intensity. */

    float cal_min_;             /**< The minimum display intensity. */

    float slice_duration_;      /**< The time to acquire 1 MRI slice. */

    float toffset_;             /**< The time axis shift. */

    int   glmax_;               /**< Unused variable. */

    int   glmin_;               /**< Unused variable. */

    char  descrip_[80];         /**< An image description text. */

    char  aux_file_[24];        /**< An auxiliary filename. */

    short qform_code_ ;         /**< The NIFTI_XFORM_* code. */

    short sform_code_ ;         /**< The NIFTI_XFORM_* code. */

    float quatern_b_;           /**< The quaternion b param to define rotation matrix. */

    float quatern_c_;           /**< The quaternion c param to define rotation matrix. */

    float quatern_d_;           /**< The quaternion d param to define rotation matrix. */

    float qoffset_x_;           /**< The quaternion x shift a.k.a offset of image. */

    float qoffset_y_;           /**< The quaternion y shift a.k.a offset of image. */

    float qoffset_z_;           /**< The quaternion z shift a.k.a offset of image. */

    float srow_x_[4];           /**< The 1st row affine transform. */

    float srow_y_[4];           /**< The 2nd row affine transform. */

    float srow_z_[4];           /**< The 3rd row affine transform. */

    char intent_name_[16];      /**< The 'name' or meaning of the data.  */

    char magic_[4];             /**< The image data location . It must be "ni1\0" for image data stored in [.img] file
                                 or "n+1\0" if image data is stored in the same file as the header [.nii]. */

} Nifti1Header;


/**
   \brief Stream output of the Nifti1Header.
   \param [out] out Stream for output.
   \param [in] header The Nifti1Header to be flushed in the output.
   \return [std::ostream] Prints the Nifti1Header to be flushed in the output stream.
*/
inline std::ostream & operator << (std::ostream &out, const Nifti1Header &header);


/**
 * \struct Nifti2Header
 * \brief Structure collecting the attributes of the 540 bytes header of a voxelized image in NifTI-2 format [.nii].
 */

typedef struct Nifti2Header{

    /**
     * \brief The Nifti2Header constructor.
     */
    inline Nifti2Header() : hdr_size_(0), magic_{'\0'}, datatype_(0), bitpix_(0), dim_{0}, intent_p1_(0.), intent_p2_(0.),
        intent_p3_(0.), pixdim_{0.}, vox_offset_(0), scl_slope_(0.), scl_inter_(0.), cal_max_(0.), cal_min_(0.),
        slice_duration_(0.), toffset_(0.), slice_start_(0), slice_end_(0), descrip_{'\0'}, aux_file_{'\0'}, qform_code_(0),
        sform_code_(0), quatern_b_(0.), quatern_c_(0.), quatern_d_(0.), qoffset_x_(0.), qoffset_y_(0.), qoffset_z_(0.),
        srow_x_{0.}, srow_y_{0.}, srow_z_{0.}, slice_code_(0), xyzt_units_(0), intent_code_(0), intent_name_{'\0'},
        dim_info_{'\0'}, unused_str_{'\0'}
    {}

    int   hdr_size_;            /**< Size of the header in bytes. It must be 540. */

    char  magic_[8];            /**< MUST be valid signature. */

    int16_t datatype_;          /**< Defines the data type of the NifTI-2 image. */

    int16_t bitpix_;            /**< The number of bits per voxel. */

    int64_t dim_[8];            /**< The image dimensions. */

    double intent_p1_;          /**< The 1st intent parameter. */

    double intent_p2_;          /**< The 2nd intent parameter. */

    double intent_p3_;          /**< The 3rd intent parameter. */

    double pixdim_[8];          /**< The image spacings. */

    int64_t vox_offset_;        /**< The offset into the .nii file where the image data start. */

    double scl_slope_;          /**< The slope of the data scaling. */

    double scl_inter_;          /**< The offset of the data scaling. */

    double cal_max_;            /**< The maximum display intensity. */

    double cal_min_;            /**< The minimum display intensity. */

    double slice_duration_;     /**< Time for 1 slice. */

    double toffset_;            /**< The time axis shift. */

    int64_t slice_start_;       /**< The index of the first slice. */

    int64_t slice_end_;         /**< The index of the last slice. */

    char  descrip_[80];         /**< An image description text. */

    char  aux_file_[24];        /**< An auxiliary filename. */

    int qform_code_;            /**< NIFTI_XFORM_* code. */

    int sform_code_;            /**< NIFTI_XFORM_* code. */

    double quatern_b_;          /**< The quaternion b param to define rotation matrix. */

    double quatern_c_;          /**< The quaternion c param to define rotation matrix. */

    double quatern_d_;          /**< The quaternion d param to define rotation matrix. */

    double qoffset_x_;          /**< The quaternion x shift a.k.a offset of image. */

    double qoffset_y_;          /**< The quaternion y shift a.k.a offset of image. */

    double qoffset_z_;          /**< The quaternion z shift a.k.a offset of image. */

    double srow_x_[4];          /**< The 1st row affine transform. */

    double srow_y_[4];          /**< The 2nd row affine transform. */

    double srow_z_[4];          /**< The 3rd row affine transform. */

    int slice_code_;            /**< Slice timing order.   */

    int xyzt_units_;            /**< The units of the image spacings -> pixdim[1..4]. */

    int intent_code_;           /**< NIFTI_INTENT_* code.  */

    char intent_name_[16];      /**< The 'name' or meaning of the data.  */

    char dim_info_;             /**< The MRI slice ordering.   */

    char unused_str_[15];       /**< unused, filled with \0 */

} Nifti2Header;


/** \} End of Doxygen Groups*/
} //end of namespace IMP

#endif //IMP_ENGINE_IMAGE_IO_NII_IO_HPP_

#include "IMP/engine/image_io/nii_io.tpp"
