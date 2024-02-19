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

#ifndef IMP_ENGINE_IMAGE_VOXEL_IMAGE_TPP_
#define IMP_ENGINE_IMAGE_VOXEL_IMAGE_TPP_

#include "IMP/engine/image/voxel_image.hpp"
#include "IMP/engine/image_io/h33_io.hpp"
#include "IMP/engine/image_io/inr_io.hpp"
#include "IMP/engine/image_io/mha_io.hpp"
#include "IMP/engine/image_io/mhd_io.hpp"
#include "IMP/engine/image_io/nii_io.hpp"



namespace IMP {

template <typename VOXELTYPE>
VoxelImage<VOXELTYPE>::VoxelImage() : data_(), dimensions_({-1, -1, -1}), spacing_({-1., -1., -1.}), origin_({0., 0., 0.}),
                                      num_voxels_(-1), memory_size_(-1), compressed_bytes_(-1),
                                      datatype_(ImageDataType::no_type), compression_state_(ImageCompression::OFF)
{}


template <typename VOXELTYPE>
VoxelImage<VOXELTYPE>::VoxelImage(const VoxelImage<VOXELTYPE> &image) : data_(image.data_),  dimensions_(image.dimensions_),
                                                                        spacing_(image.spacing_), origin_(image.origin_),
                                                                        num_voxels_(image.num_voxels_), memory_size_(image.memory_size_),
                                                                        datatype_(image.datatype_), compression_state_(image.compression_state_),
                                                                        compressed_bytes_(image.compressed_bytes_)
{
    // Construct voxelized image by copying the image's resources.
}


template <typename VOXELTYPE>
VoxelImage<VOXELTYPE>::VoxelImage(VoxelImage<VOXELTYPE> &&image) : data_(std::move(image.data_)),  dimensions_(std::move(image.dimensions_)),
                                                                   spacing_(std::move(image.spacing_)), origin_(std::move(image.origin_)),
                                                                   num_voxels_(image.num_voxels_), memory_size_(image.memory_size_),
                                                                   datatype_(image.datatype_), compression_state_(image.compression_state_),
                                                                   compressed_bytes_(image.compressed_bytes_)
{}


template <typename VOXELTYPE>
VoxelImage<VOXELTYPE>::~VoxelImage()
{}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::InitEmptyImage(const Vec<3, int> &dimensions,
                                           const Vec<3, double> &spacing,
                                           const Vec<3, double> &origin)
{
    // Set image properties.
    this->dimensions_ = dimensions;
    this->origin_ = origin;
    this->spacing_ = spacing;
    this->num_voxels_ = dimensions[0] * dimensions[1] * dimensions[2];
    this->memory_size_ = this->num_voxels_ * sizeof(VOXELTYPE);
    this->datatype_ = IMP::ImageDataType::no_type;
    this->compression_state_ = IMP::ImageCompression::OFF;
    this->compressed_bytes_ = -1;

    // Clear image.
    this->data_.clear();

    // Initialize image.
    this->data_.assign(this->num_voxels_, static_cast<VOXELTYPE>(0));

}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetDimensions(const int &x_dim, const int &y_dim, const int &z_dim)
{
    // Check if possitive dimensions were given.
    if (x_dim < 0 || y_dim < 0 || z_dim < 0) {
        throw std::invalid_argument(Logger::Error("VoxelImage dimensions have to be positive.").c_str());
    }

    // Set the dimensions.
    this->dimensions_.Set({x_dim, y_dim, z_dim});
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetSpacing(const double &x_spacing, const double &y_spacing, const double &z_spacing)
{
    // Check if negative voxel spacing was given.
    if (x_spacing < 0.000001 || y_spacing < 0.000001 || z_spacing < 0.000001) {
        throw std::invalid_argument(Logger::Error("Very small or negative voxel spacing value was given. "
                                                  "Image's voxel spacing has to be positive.").c_str());
    }

    // Set the voxel spacing.
    this->spacing_.Set({x_spacing, y_spacing, z_spacing});
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetOrigin(const double &x_origin, const double &y_origin, const double &z_origin)
{
    // Set the origin.
    this->origin_.Set({x_origin, y_origin, z_origin});
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetNumVoxels(const int &num_voxels)
{
    // Check if negative number of voxels was given.
    if (num_voxels < 0) {
        throw std::invalid_argument(Logger::Error("Cannot set number of voxels to negative number for object of VoxelImage type.").c_str());
    }

    // Set the number of voxels.
    this->num_voxels_ = num_voxels;
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetMemorySize(const int &memory_size)
{
    // Check if negative memory size was given.
    if (memory_size < 0) {
        throw std::invalid_argument(Logger::Error("Cannot set memory size to negative number for object of VoxelImage type.").c_str());
    }

    // Set the memory size.
    this->memory_size_ = memory_size;
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetDataType(const ImageDataType &datatype)
{
    // Set the datatype.
    this->datatype_ = datatype;
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetCompression(const ImageCompression &compression_state)
{
    // Set the compression state.
    this->compression_state_ = compression_state;
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SetCompressedBytes(const int &compressed_bytes)
{
    // Check if negative number of compressed bytes was given.
    if (compressed_bytes < 0) {
        throw std::invalid_argument(Logger::Error("Negative number of compressed bytes was given. "
                                                  "Image number of compressed bytes has to be positive.").c_str());
    }

    // Set the number of compressed bytes.
    this->compressed_bytes_ = compressed_bytes;
}


template <typename VOXELTYPE>
bool VoxelImage<VOXELTYPE>::IsInitialized() const
{
    // TRUE if image is initialized.
    return ((this->num_voxels_ > 0) &&
            (this->num_voxels_ == this->dimensions_[0]*this->dimensions_[1]*this->dimensions_[2]) &&
            (this->spacing_[0] > 0.000001) && (this->spacing_[1] > 0.000001) && (this->spacing_[2] > 0.000001) &&
            (this->memory_size_ == static_cast<int>(sizeof(VOXELTYPE))*this->num_voxels_) &&
            (this->datatype_ >= IMP::ImageDataType::char_type || this->datatype_ <= IMP::ImageDataType::bool_type));
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::LoadFrom(const std::string &img_filename)
{
    // // Check if there is extension in the voxelized image filename.
    // if (img_filename.find(".") == std::string::npos) {
    //     throw std::invalid_argument(Logger::Error("Cannot load VoxelImage to file with no extension. "
    //                                               "Supported format: [.h33 | .hdr | .inr | .mha | .mhd | .nii] "
    //                                               "Check given filename: " + img_filename).c_str());
    // }

    // // Get the extension of the voxelized image filename.
    // auto extension = img_filename.substr(img_filename.find_last_of("."));
    // if (extension == ".gz") {
    //     // Stripped out the compression extension.
    //     auto stripped_filename = img_filename.substr(0, img_filename.find_last_of("."));

    //     // Get the pre-extension.
    //     auto pre_ext = stripped_filename.substr(stripped_filename.find_last_of("."));

    //     // Update the extension.
    //     extension = pre_ext + extension;
    // }

    // // Load the appropriate format of the voxelized image according to it's extension.
    // if (extension == ".hdr" || extension == ".hdr.gz" || extension == ".nii" || extension == ".nii.gz") {
    //     NiiIO nii_io;
    //     nii_io.Load<VOXELTYPE>(img_filename.c_str(), *this);
    // }
    // else if (extension == ".mha") {
    //     MhaIO mha_io;
    //     mha_io.Load<VOXELTYPE>(img_filename.c_str(), *this);
    // }
    // else if (extension == ".mhd") {
    //     MhdIO mhd_io;
    //     mhd_io.Load<VOXELTYPE>(img_filename.c_str(), *this);
    // }
    // else if (extension == ".h33") {
    //     H33IO h33_io;
    //     h33_io.Load<VOXELTYPE>(img_filename.c_str(), *this);
    // }
    // else if (extension == ".inr" || extension == ".inr.gz") {
    //     InrIO inr_io;
    //     inr_io.Load<VOXELTYPE>(img_filename.c_str(), *this);
    // }
    // else {
    //     throw std::invalid_argument(Logger::Error("Cannot load VoxelImage. "
    //                                               "Supported format: [.h33 | .hdr | .inr | .mha | .mhd | .nii] "
    //                                               "Check given filename: " + img_filename).c_str());
    // }
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::SaveTo(const std::string &img_filename, ImageCompression compression)
{
    // // Check if there is extension in the voxelized image filename.
    // if (img_filename.find(".") == std::string::npos) {
    //     throw std::invalid_argument(Logger::Error("Cannot save VoxelImage to file with no extension. "
    //                                               "Supported format: [.h33 | .hdr | .inr | .mha | .mhd | .nii] "
    //                                               "Check given filename: " + img_filename).c_str());
    // }

    // // Get the extension of the voxelized image filename.
    // auto extension = img_filename.substr(img_filename.find_last_of("."));
    // if (extension == ".gz") {
    //     // Stripped out the compression extension.
    //     auto stripped_filename = img_filename.substr(0, img_filename.find_last_of("."));

    //     // Get the pre-extension.
    //     auto pre_ext = stripped_filename.substr(stripped_filename.find_last_of("."));

    //     // Update the extension.
    //     extension = pre_ext + extension;
    // }

    // // Save to corresponding format regarding the extension of the given image filename.
    // if (extension == ".hdr" || extension == ".hdr.gz" || extension == ".nii" || extension == ".nii.gz") {
    //     NiiIO nii_io;
    //     nii_io.Save<VOXELTYPE>(*this, img_filename.c_str(), compression);
    // }
    // else if (extension == ".mha") {
    //     MhaIO mha_io;
    //     mha_io.Save<VOXELTYPE>(*this, img_filename.c_str(), compression);
    // }
    // else if (extension == ".mhd") {
    //     MhdIO mhd_io;
    //     mhd_io.Save<VOXELTYPE>(*this, img_filename.c_str(), compression);
    // }
    // else if (extension == ".h33") {
    //     H33IO h33_io;
    //     h33_io.Save<VOXELTYPE>(*this, img_filename.c_str(), compression);
    // }
    // else if (extension == ".inr" || extension == ".inr.gz") {
    //     InrIO inr_io;
    //     inr_io.Save<VOXELTYPE>(*this, img_filename.c_str(), compression);
    // }
    // else {
    //     throw std::invalid_argument(Logger::Error("Cannot load VoxelImage. "
    //                                               "Supported format: [.h33 | .hdr | .inr | .mha | .mhd | .nii] "
    //                                               "Check given filename: " + img_filename).c_str());
    // }

}


template <typename VOXELTYPE>
bool VoxelImage<VOXELTYPE>::operator == (const VoxelImage<VOXELTYPE> &image) const
{
    // Compare voxelized images for equality.
    return ((this->data_ == image.data_) &&
            (this->dimensions_ == image.dimensions_) &&
            (this->spacing_ == image.spacing_) &&
            (this->origin_ == image.origin_) &&
            (this->num_voxels_ == image.num_voxels_) &&
            (this->memory_size_ == image.memory_size_) &&
            (this->datatype_ == image.datatype_) &&
            (this->compression_state_ == image.compression_state_) &&
            (this->compressed_bytes_ == image.compressed_bytes_)
           );
}


template <typename VOXELTYPE>
bool VoxelImage<VOXELTYPE>::operator != (const VoxelImage<VOXELTYPE> &image) const
{
    // Compare voxelized images for inequality.
    return !(*this == image);
}


template <typename VOXELTYPE>
VoxelImage<VOXELTYPE> & VoxelImage<VOXELTYPE>::operator = (const VoxelImage<VOXELTYPE> &image)
{
    if (this != &image) {
        // Assign by copying image resources.
        this->data_ = image.data_;
        this->dimensions_ = image.dimensions_;
        this->spacing_ = image.spacing_;
        this->origin_ = image.origin_;
        this->num_voxels_ = image.num_voxels_;
        this->memory_size_ = image.memory_size_;
        this->datatype_ = image.datatype_;
        this->compression_state_ == image.compression_state_;
        this->compressed_bytes_ == image.compressed_bytes_;
    }
    return *this;
}


template <typename VOXELTYPE>
VoxelImage<VOXELTYPE> & VoxelImage<VOXELTYPE>::operator = (VoxelImage<VOXELTYPE> &&image)
{
    if (this != &image) {
        // Assign by moving image resources.
        this->data_ = std::move(image.data_);
        this->dimensions_ = std::move(image.dimensions_);
        this->spacing_ = std::move(image.spacing_);
        this->origin_ = std::move(image.origin_);
        this->num_voxels_ = std::move(image.num_voxels_);
        this->memory_size_ = std::move(image.memory_size_);
        this->datatype_ = std::move(image.datatype_);
        this->compression_state_ == std::move(image.compression_state_);
        this->compressed_bytes_ == std::move(image.compressed_bytes_);

        // Reset image resources.
        image.data_.clear();
        image.dimensions_.Set({0, 0, 0});
        image.dimensions_.Set({-1, -1, -1});
        image.spacing_.Set({-1., -1., -1.});
        image.origin_.Set({0., 0., 0.});
        image.num_voxels_ = -1;
        image.memory_size_ = -1;
        image.datatype_ = ImageDataType::no_type;
        image.compression_state_ = ImageCompression::OFF;
        image.compressed_bytes_ = -1;
    }
    return *this;
}


template <typename VOXELTYPE>
VOXELTYPE VoxelImage<VOXELTYPE>::operator [] (const int &id) const
{
    return this->data_.at(id);
}


template <typename VOXELTYPE>
VOXELTYPE VoxelImage<VOXELTYPE>::operator () (const int &i, const int &j, const int &k) const
{
    int id = i + (j*this->dimensions_[0]) + (k*this->dimensions_[0]*this->dimensions_[1]);

    return this->data_.at(id);
}


template <typename VOXELTYPE>
void VoxelImage<VOXELTYPE>::Clear()
{
    this->dimensions_.Set({-1, -1, -1});
    this->origin_.Set({0., 0., 0.});
    this>spacing_.Set({-1.,-1.,-1.});
    this->num_voxels_ = -1;
    this->memory_size_ = -1;
    this->datatype_ = ImageDataType::no_type;
    this->compression_state_ = ImageCompression::OFF;
    this->compressed_bytes_ = -1;
    this->data_.clear();
}


} // End of namespace IMP

#endif //IMP_ENGINE_IMAGE_VOXEL_IMAGE_TPP_
