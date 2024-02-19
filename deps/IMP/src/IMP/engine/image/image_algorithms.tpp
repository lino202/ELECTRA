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


#ifndef IMP_ENGINE_IMAGE_IMAGE_ALGORITHMS_TPP_
#define IMP_ENGINE_IMAGE_IMAGE_ALGORITHMS_TPP_


#include "IMP/engine/image/image_algorithms.hpp"


namespace IMP {


// template<typename VOXELTYPE>
// CGAL::Image_3 CgalImage_3FromIMPVoxelImage(const VoxelImage<VOXELTYPE> &vox_image)
// {
//     // Check if voxelized image is initialized before conversion.
//     if (!vox_image.IsInitialized()) {
//         throw std::runtime_error(Logger::Error("Cannot convert IMP::VoxelImage to CGAL::Image_3."
//                                                " IMP::VoxelImage is not initialized.").c_str());
//     }

//     // Initialize point_image pointer to copy the voxelized image data.
//     _image *image = ::_initImage();

//     // Set image dimensions.
//     image->xdim = vox_image.Dimensions()[0];
//     image->ydim = vox_image.Dimensions()[1];
//     image->zdim = vox_image.Dimensions()[2];

//     // Set image voxel spacing.
//     image->vx = vox_image.Spacing()[0];
//     image->vy = vox_image.Spacing()[1];
//     image->vz = vox_image.Spacing()[2];

//     image->cx = vox_image.Origin()[0];
//     image->cy = vox_image.Origin()[1];
//     image->cz = vox_image.Origin()[2];

//     // Set image value type to scalar.
//     image->vectMode = VM_SCALAR;

//     // Set image vectorial dimension.
//     image->vdim = 1;

//     // Set image endianness.
//     image->endianness = ::_getEndianness();

//     // Set image word dimension.
//     image->wdim = sizeof(VOXELTYPE);

//     // Set image datatype and sign.
//     if (vox_image.DataType() == ImageDataType::bool_type) {
//         image->wordKind = WK_FIXED;
//         image->sign = SGN_UNKNOWN;
//     }
//     else if (vox_image.DataType() == ImageDataType::char_type ||
//              vox_image.DataType() == ImageDataType::short_int_type ||
//              vox_image.DataType() == ImageDataType::int_type ||
//              vox_image.DataType() == ImageDataType::long_int_type) {
//         image->wordKind = WK_FIXED;
//         image->sign = SGN_SIGNED;
//     }
//     else if (vox_image.DataType() == ImageDataType::uchar_type ||
//              vox_image.DataType() == ImageDataType::ushort_int_type ||
//              vox_image.DataType() == ImageDataType::uint_type ||
//              vox_image.DataType() == ImageDataType::ulong_int_type) {
//         image->wordKind = WK_FIXED;
//         image->sign = SGN_UNSIGNED;
//     }
//     else if (vox_image.DataType() == ImageDataType::float_type ||
//              vox_image.DataType() == ImageDataType::double_type){
//         image->wordKind = WK_FLOAT;
//         image->sign = SGN_UNKNOWN;
//     }

//     // Allocate and copy image data.
//     image->data = ::ImageIO_alloc(vox_image.MemorySize());
//     std::memcpy(image->data, static_cast<void*>(&std::vector<VOXELTYPE>(std::move(vox_image.Data()))[0]), vox_image.MemorySize());

//     // Return the cgal image.
//     return CGAL::Image_3(image);

// }





} // End of namespace IMP

#endif //IMP_ENGINE_IMAGE_IMAGE_ALGORITHMS_TPP_
