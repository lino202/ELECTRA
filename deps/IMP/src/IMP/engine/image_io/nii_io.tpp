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


#ifndef IMP_ENGINE_IMAGE_IO_NII_IO_TPP_
#define IMP_ENGINE_IMAGE_IO_NII_IO_TPP_


#include "IMP/engine/image_io/nii_io.hpp"


namespace IMP {


NiiIO::NiiIO()
{}


NiiIO::~NiiIO()
{}


template<typename VOXELTYPE>
void NiiIO::Load(const std::string &filename, VoxelImage<VOXELTYPE> &image)
{
    // The NifTI image string stream.
    std::stringstream nii_sstream(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

    // Open NifTI image according to the extension.
    std::string extension = filename.substr(filename.find_last_of("."));
    if (extension == ".hdr" || extension == ".nii") {
        // Open NifTI file.
        std::ifstream nii_file(filename.c_str(), std::ios_base::in | std::ios_base::binary);

        // Check if file is open.
        if (!nii_file.is_open()) {
            throw std::runtime_error(Logger::Error("Cannot open as NifTI [.hdr | .nii]. "
                                                   "Check file: " + filename).c_str());
        }

        // Read the NifTI file data.
        nii_sstream << nii_file.rdbuf();
        nii_file.close();
    }
    else {
        // Search for compressed file extension.
        extension = filename.substr(filename.find_last_of(".")-4);
        if (extension == ".hdr.gz" || extension == ".nii.gz") {
            // Open compressed NifTI file.
            std::ifstream nii_comp_file(filename.c_str(), std::ios_base::in | std::ios_base::binary);

            // Check if file is open.
            if (!nii_comp_file.is_open()) {
                throw std::runtime_error(Logger::Error("Cannot open as compressed NifTI [.hdr.gz | .nii.gz]. "
                                                       "Check file: " + filename).c_str());
            }

            // Header decompressing scope.
            {
                // Read compressed binary data.
                boost::iostreams::filtering_istream input;
                input.push(boost::iostreams::gzip_decompressor());
                input.push(nii_comp_file);
                boost::iostreams::copy(input, nii_sstream);
            } // end of Header decompressing scope.

            // Close compressed NifTI file.
            nii_comp_file.close();
        }
        else {
            throw std::invalid_argument(Logger::Error("Cannot load file of unknown format as NifTI "
                                                      "[.hdr | .hdr.gz | .nii | .nii.gz]. Check file: " + filename).c_str());
        }
    }


    // Read the header size value.
    int header_size = 0;
    nii_sstream.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));

    // Choose to load NifTI-1 or NifTI-2 according to header size.
    switch (header_size) {
    case 348:
    {
        // Read header in NifTI-1 format.
        Nifti1Header ni1_header;

        // Assign the header size.
        ni1_header.hdr_size_ = header_size;

        // Read the rest header attributes.
        nii_sstream.read(ni1_header.data_type_, sizeof(ni1_header.data_type_));
        nii_sstream.read(ni1_header.db_name_, sizeof(ni1_header.db_name_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.extents_), sizeof(ni1_header.extents_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.session_error_), sizeof(ni1_header.session_error_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.regular_), sizeof(ni1_header.regular_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.dim_info_), sizeof(ni1_header.dim_info_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.dim_), sizeof(ni1_header.dim_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.intent_p1_), sizeof(ni1_header.intent_p1_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.intent_p2_), sizeof(ni1_header.intent_p2_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.intent_p3_), sizeof(ni1_header.intent_p3_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.intent_code_), sizeof(ni1_header.intent_code_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.datatype_), sizeof(ni1_header.datatype_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.bitpix_), sizeof(ni1_header.bitpix_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.slice_start_), sizeof(ni1_header.slice_start_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.pixdim_), sizeof(ni1_header.pixdim_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.vox_offset_), sizeof(ni1_header.vox_offset_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.scl_slope_), sizeof(ni1_header.scl_slope_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.scl_inter_), sizeof(ni1_header.scl_inter_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.slice_end_), sizeof(ni1_header.slice_end_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.slice_code_), sizeof(ni1_header.slice_code_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.xyzt_units_), sizeof(ni1_header.xyzt_units_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.cal_max_), sizeof(ni1_header.cal_max_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.cal_min_), sizeof(ni1_header.cal_min_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.slice_duration_), sizeof(ni1_header.slice_duration_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.toffset_), sizeof(ni1_header.toffset_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.glmax_), sizeof(ni1_header.glmax_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.glmin_), sizeof(ni1_header.glmin_));
        nii_sstream.read(ni1_header.descrip_, sizeof(ni1_header.descrip_));
        nii_sstream.read(ni1_header.aux_file_, sizeof(ni1_header.aux_file_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.qform_code_), sizeof(ni1_header.qform_code_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.sform_code_), sizeof(ni1_header.sform_code_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.quatern_b_), sizeof(ni1_header.quatern_b_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.quatern_c_), sizeof(ni1_header.quatern_c_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.quatern_d_), sizeof(ni1_header.quatern_d_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.qoffset_x_), sizeof(ni1_header.qoffset_x_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.qoffset_y_), sizeof(ni1_header.qoffset_y_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.qoffset_z_), sizeof(ni1_header.qoffset_z_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.srow_x_), sizeof(ni1_header.srow_x_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.srow_y_), sizeof(ni1_header.srow_y_));
        nii_sstream.read(reinterpret_cast<char*>(&ni1_header.srow_z_), sizeof(ni1_header.srow_z_));
        nii_sstream.read(ni1_header.intent_name_, sizeof(ni1_header.intent_name_));
        nii_sstream.read(ni1_header.magic_, sizeof(ni1_header.magic_));

        // Check for image dimensions.
        if (ni1_header.dim_[0] != 3) {
            throw std::runtime_error(Logger::Error("Cannot load NifTI-1 image if has not 3 dimensions. "
                                                   "Check file: " + filename).c_str());
        }

        // Check magic number
        if (ni1_header.magic_[0] != 'n' &&
            ni1_header.magic_[1] != 'i' && ni1_header.magic_[1] != '+' &&
            ni1_header.magic_[2] != '1' && ni1_header.magic_[3] != '\0') {
            throw std::runtime_error(Logger::Error("Cannot load NifTI-1 image of unknown magic number."
                                                   " Check file: " + filename).c_str());
        }

        // Initialize image.
        image.SetDimensions(static_cast<int>(ni1_header.dim_[1]),
                            static_cast<int>(ni1_header.dim_[2]),
                            static_cast<int>(ni1_header.dim_[3]));

        image.SetSpacing(static_cast<double>(ni1_header.pixdim_[1]),
                         static_cast<double>(ni1_header.pixdim_[2]),
                         static_cast<double>(ni1_header.pixdim_[3]));

        // Origin correction for coordinates system.
        image.SetOrigin(static_cast<double>(-1.f*ni1_header.qoffset_x_),
                        static_cast<double>(-1.f*ni1_header.qoffset_y_),
                        static_cast<double>(ni1_header.qoffset_z_));

        image.SetNumVoxels(static_cast<int>(ni1_header.dim_[1]*ni1_header.dim_[2]*ni1_header.dim_[3]));

        image.SetMemorySize(static_cast<int>(sizeof(VOXELTYPE))*image.NumVoxels());

        // Initialize the data type of the image.
        if (ni1_header.datatype_ == 256) { image.SetDataType(ImageDataType::char_type); }
        else if (ni1_header.datatype_ == 2) { image.SetDataType(ImageDataType::uchar_type); }
        else if (ni1_header.datatype_ == 4) { image.SetDataType(ImageDataType::short_int_type); }
        else if (ni1_header.datatype_ == 512) { image.SetDataType(ImageDataType::ushort_int_type); }
        else if (ni1_header.datatype_ == 8) { image.SetDataType(ImageDataType::int_type); }
        else if (ni1_header.datatype_ == 768) { image.SetDataType(ImageDataType::uint_type); }
        else if (ni1_header.datatype_ == 1024) { image.SetDataType(ImageDataType::long_int_type); }
        else if (ni1_header.datatype_ == 1280) { image.SetDataType(ImageDataType::ulong_int_type); }
        else if (ni1_header.datatype_ == 16) { image.SetDataType(ImageDataType::float_type); }
        else if (ni1_header.datatype_ == 64) { image.SetDataType(ImageDataType::double_type); }
        else { throw std::runtime_error(Logger::Error("Cannot load NifTI-1 image of unsupported data type."
                                                      " Check file: " + filename).c_str()); }

        // Check image data type compatibility.
        if ( !(image.DataType() == ImageDataType::char_type && sizeof(VOXELTYPE) == sizeof(char)) &&
             !(image.DataType() == ImageDataType::uchar_type && sizeof(VOXELTYPE) == sizeof(unsigned char)) &&
             !(image.DataType() == ImageDataType::short_int_type && sizeof(VOXELTYPE) == sizeof(short)) &&
             !(image.DataType() == ImageDataType::ushort_int_type && sizeof(VOXELTYPE) == sizeof(unsigned short)) &&
             !(image.DataType() == ImageDataType::int_type && sizeof(VOXELTYPE) == sizeof(int)) &&
             !(image.DataType() == ImageDataType::uint_type && sizeof(VOXELTYPE) == sizeof(unsigned int)) &&
             !(image.DataType() == ImageDataType::long_int_type && sizeof(VOXELTYPE) == sizeof(long int)) &&
             !(image.DataType() == ImageDataType::ulong_int_type && sizeof(VOXELTYPE) == sizeof(unsigned long int)) &&
             !(image.DataType() == ImageDataType::float_type && sizeof(VOXELTYPE) == sizeof(float)) &&
             !(image.DataType() == ImageDataType::double_type && sizeof(VOXELTYPE) == sizeof(double))
             ) { throw std::runtime_error(Logger::Error("Cannot load NifTI-1 image in VoxelImage of different data type."
                                                        " Check file: " + filename).c_str()); }

        // Reset image data container.
        image.EditData().clear();
        image.EditData().assign(image.NumVoxels(), 0);

        // Read the image data.
        if (extension == ".hdr") {
            // Open the binary image (.img) file.
            std::string img_filename = filename.substr(0, filename.find_last_of(".")) + ".img";
            std::ifstream img_file(img_filename.c_str(), std::ios_base::in | std::ios_base::binary);

            // Check if file is open.
            if (!img_file.is_open()) {
                throw std::runtime_error(Logger::Error("Cannot open the NifTI-1 image file [.img]."
                                                       " Check file: " + img_filename).c_str());
            }

            // Create binary data buffer.
            VOXELTYPE *buffer = new VOXELTYPE [image.NumVoxels()];

            // Read binary image data.
            img_file.read((char *)buffer, image.MemorySize());

            // Store the binary data in the voxelized image's data vector.
            for(auto i = 0; i != image.NumVoxels(); ++i) {
                image.EditData()[i] = static_cast<VOXELTYPE>(buffer[i]);
            }

            // Clean-up buffer allocated memory.
            delete[] buffer;
            buffer = nullptr;
        }
        else if (extension == ".hdr.gz") {
            // Open the compressed binary image (.img) file.
            std::string img_comp_filename = filename.substr(0, filename.find_last_of(".")-4) + ".img.gz";
            std::ifstream img_comp_file(img_comp_filename.c_str(), std::ios_base::in | std::ios_base::binary);

            // Check if file is open.
            if (!img_comp_file.is_open()) {
                throw std::runtime_error(Logger::Error("Cannot open the compressed NifTI-1 image file [.img.gz]."
                                                       " Check file: " + img_comp_filename).c_str());
            }

            // Compressed file reading scope.
            {
                // Read compressed binary data.
                boost::iostreams::filtering_istream input;
                input.push(boost::iostreams::gzip_decompressor());
                input.push(img_comp_file);

                // Get character pointer to the image data vector.
                char *data_ptr = reinterpret_cast<char*>(image.EditData().data());

                // Assign data to the image data vector byte by byte.
                for (int i = 0; i != image.MemorySize(); ++i) {
                    input.read(&data_ptr[i], 1);
                }

                // Release character pointer.
                data_ptr = nullptr;

            } // End of compressed file reading scope.

        }
        else if (extension == ".nii" || extension == ".nii.gz") {

            // Read the image data starting at the vox_offset byte of the file.
            nii_sstream.seekg(ni1_header.vox_offset_, std::ios::beg);

            // Create binary data buffer.
            VOXELTYPE *buffer = new VOXELTYPE [image.NumVoxels()];

            // Read binary image data.
            nii_sstream.read((char *)buffer, image.MemorySize());

            // Store the binary data in the voxelized image's data vector.
            for(auto i = 0; i != image.NumVoxels(); ++i) {
                image.EditData()[i] = static_cast<VOXELTYPE>(buffer[i]);
            }

            // Clean-up buffer allocated memory.
            delete[] buffer;
            buffer = nullptr;
        }

        // Break case NifTI-1 format.
        break;
    }
    case 540:
    {
        // Read header in NifTI-2 format.
        Nifti2Header ni2_header;

        // Assign the header size.
        ni2_header.hdr_size_ = header_size;

        // Read the rest header attributes.
        nii_sstream.read(ni2_header.magic_, sizeof(ni2_header.magic_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.datatype_), sizeof(ni2_header.datatype_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.bitpix_), sizeof(ni2_header.bitpix_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.dim_), sizeof(ni2_header.dim_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.intent_p1_), sizeof(ni2_header.intent_p1_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.intent_p2_), sizeof(ni2_header.intent_p2_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.intent_p3_), sizeof(ni2_header.intent_p3_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.pixdim_), sizeof(ni2_header.pixdim_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.vox_offset_), sizeof(ni2_header.vox_offset_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.scl_slope_), sizeof(ni2_header.scl_slope_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.scl_inter_), sizeof(ni2_header.scl_inter_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.cal_max_), sizeof(ni2_header.cal_max_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.cal_min_), sizeof(ni2_header.cal_min_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.slice_duration_), sizeof(ni2_header.slice_duration_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.toffset_), sizeof(ni2_header.toffset_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.slice_start_), sizeof(ni2_header.slice_start_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.slice_end_), sizeof(ni2_header.slice_end_));
        nii_sstream.read(ni2_header.descrip_, sizeof(ni2_header.descrip_));
        nii_sstream.read(ni2_header.aux_file_, sizeof(ni2_header.aux_file_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.qform_code_), sizeof(ni2_header.qform_code_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.sform_code_), sizeof(ni2_header.sform_code_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.quatern_b_), sizeof(ni2_header.quatern_b_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.quatern_c_), sizeof(ni2_header.quatern_c_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.quatern_d_), sizeof(ni2_header.quatern_d_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.qoffset_x_), sizeof(ni2_header.qoffset_x_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.qoffset_y_), sizeof(ni2_header.qoffset_y_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.qoffset_z_), sizeof(ni2_header.qoffset_z_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.srow_x_), sizeof(ni2_header.srow_x_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.srow_y_), sizeof(ni2_header.srow_y_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.srow_z_), sizeof(ni2_header.srow_z_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.slice_code_), sizeof(ni2_header.slice_code_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.xyzt_units_), sizeof(ni2_header.xyzt_units_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.intent_code_), sizeof(ni2_header.intent_code_));
        nii_sstream.read(ni2_header.intent_name_, sizeof(ni2_header.intent_name_));
        nii_sstream.read(reinterpret_cast<char*>(&ni2_header.dim_info_), sizeof(ni2_header.dim_info_));
        nii_sstream.read(ni2_header.unused_str_, sizeof(ni2_header.unused_str_));

        // Check for image dimensions.
        if (ni2_header.dim_[0] > 3) {
            throw std::runtime_error(Logger::Error("Cannot load NifTI-2 image [.nii] if has not 3 dimensions."
                                                   " Check file: " + filename).c_str());
        }

        // Initialize image.
        image.SetDimensions(static_cast<int>(ni2_header.dim_[1]),
                            static_cast<int>(ni2_header.dim_[2]),
                            static_cast<int>(ni2_header.dim_[3]));

        image.SetSpacing(ni2_header.pixdim_[1], ni2_header.pixdim_[2], ni2_header.pixdim_[3]);

        // Origin correction for coordinates system.
        image.SetOrigin(-1.*ni2_header.qoffset_x_, -1.*ni2_header.qoffset_y_, ni2_header.qoffset_z_);

        image.SetNumVoxels(static_cast<int>(ni2_header.dim_[1]*ni2_header.dim_[2]*ni2_header.dim_[3]));

        image.SetMemorySize(static_cast<int>(sizeof(VOXELTYPE))*image.NumVoxels());

        // Initialize the data type of the image.
        if (ni2_header.datatype_ == 256) { image.SetDataType(ImageDataType::char_type); }
        else if (ni2_header.datatype_ == 2) { image.SetDataType(ImageDataType::uchar_type); }
        else if (ni2_header.datatype_ == 4) { image.SetDataType(ImageDataType::short_int_type); }
        else if (ni2_header.datatype_ == 512) { image.SetDataType(ImageDataType::ushort_int_type); }
        else if (ni2_header.datatype_ == 8) { image.SetDataType(ImageDataType::int_type); }
        else if (ni2_header.datatype_ == 768) { image.SetDataType(ImageDataType::uint_type); }
        else if (ni2_header.datatype_ == 1024) { image.SetDataType(ImageDataType::long_int_type); }
        else if (ni2_header.datatype_ == 1280) { image.SetDataType(ImageDataType::ulong_int_type); }
        else if (ni2_header.datatype_ == 16) { image.SetDataType(ImageDataType::float_type); }
        else if (ni2_header.datatype_ == 64) { image.SetDataType(ImageDataType::double_type); }
        else { throw std::runtime_error(Logger::Error("Cannot load NifTI-2 image of unsupported data type."
                                                      " Check file: " + filename).c_str()); }

        // Check image data type compatibility.
        if ( !(image.DataType() == ImageDataType::char_type && sizeof(VOXELTYPE) == sizeof(char)) &&
             !(image.DataType() == ImageDataType::uchar_type && sizeof(VOXELTYPE) == sizeof(unsigned char)) &&
             !(image.DataType() == ImageDataType::short_int_type && sizeof(VOXELTYPE) == sizeof(short)) &&
             !(image.DataType() == ImageDataType::ushort_int_type && sizeof(VOXELTYPE) == sizeof(unsigned short)) &&
             !(image.DataType() == ImageDataType::int_type && sizeof(VOXELTYPE) == sizeof(int)) &&
             !(image.DataType() == ImageDataType::uint_type && sizeof(VOXELTYPE) == sizeof(unsigned int)) &&
             !(image.DataType() == ImageDataType::long_int_type && sizeof(VOXELTYPE) == sizeof(long int)) &&
             !(image.DataType() == ImageDataType::ulong_int_type && sizeof(VOXELTYPE) == sizeof(unsigned long int)) &&
             !(image.DataType() == ImageDataType::float_type && sizeof(VOXELTYPE) == sizeof(float)) &&
             !(image.DataType() == ImageDataType::double_type && sizeof(VOXELTYPE) == sizeof(double))
             ) { throw std::runtime_error(Logger::Error("Cannot load NifTI-2 image in VoxelImage of different data type."
                                                        " Check file: " + filename).c_str()); }

        // Reset image data container.
        image.EditData().clear();
        image.EditData().assign(image.NumVoxels(), 0);

        // Read the image data starting at the vox_offset byte of the file.
        nii_sstream.seekg(ni2_header.vox_offset_, std::ios::beg);

        // Create binary data buffer.
        VOXELTYPE *buffer = new VOXELTYPE [image.NumVoxels()];

        // Read binary image data.
        nii_sstream.read((char *)buffer, image.MemorySize());

        // Store the binary data in the voxelized image's data vector.
        for(auto i = 0; i != image.NumVoxels(); ++i) {
            image.EditData()[i] = static_cast<VOXELTYPE>(buffer[i]);
        }

        // Clean-up buffer allocated memory.
        delete[] buffer;
        buffer = nullptr;

        // Break case NifTI-2 format.
        break;
    }
    default:
        throw std::invalid_argument(Logger::Error("The header size: " + std::to_string(header_size) +
                                                  " does not comply with NifTI-1 [348] or NifTI-2 [540] "
                                                  "expected header size. Check file: " + filename).c_str());

    } // end of switch to read NifTI-1 or Nifti-2 file.

}


template<typename VOXELTYPE>
void NiiIO::Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression)
{
    // Initialize the path of the output file.
    std::string path = "";

    // Position of the last slash in the output file's name.
    std::size_t last_slash = filename.find_last_of("/\\");

    // Get the path directory of the output file.
    if (last_slash != std::string::npos) { path = filename.substr(0, last_slash); }

    // Create the path's directory if it doesn't exist.
    boost::filesystem::path dir(path);
    if (!path.empty() && !boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
    }

    // Get NifTI-1 output file extension.
    std::string extension = filename.substr(filename.find_last_of("."));
    if (extension != ".nii" && extension != ".hdr") {
        // Search for compressed file extension.
        extension = filename.substr(filename.find_last_of(".")-4);
        if (extension != ".nii.gz" && extension != ".hdr.gz") {
            throw std::invalid_argument(Logger::Error("Cannot save NifTI-1 image to unsupported format. "
                                                      "Supported format: [.nii | .nii.gz | .hdr | .hdr.gz]."
                                                      " Check output filename: " + filename).c_str());
        }
    }

    // Initialize the NifTI-1 header.
    Nifti1Header ni1_header;

    // Write the NifTI-1 header.
    ni1_header.hdr_size_ = 348;
    ni1_header.regular_ = 'r';
    ni1_header.dim_[0] = 3;
    ni1_header.dim_[1] = image.Dimensions()[0];
    ni1_header.dim_[2] = image.Dimensions()[1];
    ni1_header.dim_[3] = image.Dimensions()[2];
    ni1_header.dim_[4] = 1;
    ni1_header.dim_[5] = 1;
    ni1_header.dim_[6] = 1;
    ni1_header.dim_[7] = 1;

    if (image.DataType() == ImageDataType::char_type) { ni1_header.datatype_ = 256; ni1_header.bitpix_ = 8*sizeof(char); }
    else if (image.DataType() == ImageDataType::uchar_type) { ni1_header.datatype_ = 2; ni1_header.bitpix_ = 8*sizeof(unsigned char); }
    else if (image.DataType() == ImageDataType::short_int_type) { ni1_header.datatype_ = 4; ni1_header.bitpix_ = 8*sizeof(short int); }
    else if (image.DataType() == ImageDataType::ushort_int_type) { ni1_header.datatype_ = 512; ni1_header.bitpix_ = 8*sizeof(unsigned short int); }
    else if (image.DataType() == ImageDataType::int_type) { ni1_header.datatype_ = 8; ni1_header.bitpix_ = 8*sizeof(int); }
    else if (image.DataType() == ImageDataType::uint_type) { ni1_header.datatype_ = 768; ni1_header.bitpix_ = 8*sizeof(unsigned int); }
    else if (image.DataType() == ImageDataType::long_int_type) { ni1_header.datatype_ = 1024; ni1_header.bitpix_ = 8*sizeof(long int); }
    else if (image.DataType() == ImageDataType::ulong_int_type) { ni1_header.datatype_ = 1280; ni1_header.bitpix_ = 8*sizeof(unsigned long int); }
    else if (image.DataType() == ImageDataType::float_type) { ni1_header.datatype_ = 16; ni1_header.bitpix_ = 8*sizeof(float); }
    else if (image.DataType() == ImageDataType::double_type) { ni1_header.datatype_ = 64; ni1_header.bitpix_ = 8*sizeof(double); }
    else if (image.DataType() == ImageDataType::bool_type) { ni1_header.datatype_ = 2; ni1_header.bitpix_ = 8*sizeof(unsigned char); }
    else { throw std::runtime_error(Logger::Error("Cannot save VoxelImage of unsupported data type to NifTI-1 format.").c_str()); }

    ni1_header.pixdim_[0] = 1.f;
    ni1_header.pixdim_[1] = static_cast<float>(image.Spacing()[0]);
    ni1_header.pixdim_[2] = static_cast<float>(image.Spacing()[1]);
    ni1_header.pixdim_[3] = static_cast<float>(image.Spacing()[2]);

    if (extension == ".nii" || extension == ".nii.gz") { ni1_header.vox_offset_ = 352.f; }

    ni1_header.scl_slope_ = 1.f;
    ni1_header.qform_code_ = 2;
    ni1_header.sform_code_ = 1;
    ni1_header.quatern_d_ = 1.f;

    // Flip X,Y axis Origin of VoxelImage for NifTI-1 compatibility.
    ni1_header.qoffset_x_ = static_cast<float>(-image.Origin()[0]);
    ni1_header.qoffset_y_ = static_cast<float>(-image.Origin()[1]);
    ni1_header.qoffset_z_ = static_cast<float>(image.Origin()[2]);

    // Set affine transformation matrix elements with X, Y axis flipping.
    ni1_header.srow_x_[0] = static_cast<float>(-image.Spacing()[0]);
    ni1_header.srow_x_[3] = static_cast<float>(-image.Origin()[0]);
    ni1_header.srow_y_[1] = static_cast<float>(-image.Spacing()[1]);
    ni1_header.srow_y_[3] = static_cast<float>(-image.Origin()[1]);
    ni1_header.srow_z_[2] = static_cast<float>(image.Spacing()[2]);
    ni1_header.srow_z_[3] = static_cast<float>(image.Origin()[2]);

    // Set the NifTI-1 magic line.
    ni1_header.magic_[0] = 'n';
    if (extension == ".nii" || extension == ".nii.gz") { ni1_header.magic_[1] = '+'; }
    else { ni1_header.magic_[1] = 'i'; }
    ni1_header.magic_[2] = '1';
    ni1_header.magic_[3] = '\0';

    // Write output.
    if (extension.find("nii") != std::string::npos) { // Write in single file [.nii | .nii.gz]

        // Output stream to write file.
        std::ofstream nii_output;

        // Add ".gz" extension if missing from output file name and compressed mode is on.
        if (compression == ImageCompression::ON && extension == ".nii") {
            nii_output.open(filename + ".gz", std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }
        else {
            nii_output.open(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }

        // Check if NifTI-1 output stream file is open.
        if (!nii_output.is_open()) {
            throw std::runtime_error(Logger::Error("Could not open file to save voxelized image in "
                                                   "NifTI-1 image format [.nii | .nii.gz]. "
                                                   "Check output filename: " + filename).c_str());
        }

        // Sting stream to store uncompressed image binary data.
        std::stringstream uncompressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

        // Write the header data.
        uncompressed_data.write(reinterpret_cast<char *>(&ni1_header), sizeof(ni1_header));

        // Add four null characters padding before writing the image data
        char dummy_padding[4] = {'\0'};
        uncompressed_data.write(dummy_padding, sizeof(dummy_padding));
        uncompressed_data.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

        // Compress data if compress mode is on or [.nii.gz] extension is given.
        if (compression == ImageCompression::ON || extension == ".nii.gz") {
            // Sting stream to store the compressed image binary data.
            std::stringstream compressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

            // Compressed file writing scope.
            {
                // Create compressed output buffer.
                boost::iostreams::filtering_streambuf<boost::iostreams::input> comp_buffer;
                comp_buffer.push(boost::iostreams::gzip_compressor());
                comp_buffer.push(uncompressed_data);

                // Copy compressed output buffer in compressed binary data string stream.
                boost::iostreams::copy(comp_buffer, compressed_data);

                // Flush the data from the compressor.
                boost::iostreams::close(comp_buffer);

                // Write the compressed image binary data.
                nii_output.write(compressed_data.str().c_str(), compressed_data.str().length());

            } // End of Compressed file writing scope.
        }
        else {
            // Write the uncompressed image binary data.
            nii_output.write(uncompressed_data.str().c_str(), uncompressed_data.str().length());
        }

        // Close the NifTI-1 output stream.
        nii_output.close();

    }
    else { // Write in [.hdr | .hdr.gz] and [.img | .img.gz] files

        // Output streams to write hdr and img files.
        std::ofstream hdr_output; std::ofstream img_output;

        // Output filename of image data file.
        std::string img_filename = "";

        // Add ".gz" extension if missing from output file name and compressed mode is on.
        if (compression == ImageCompression::ON && extension == ".hdr") {
            hdr_output.open(filename + ".gz", std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

            // Change [.hdr] to [.img.gz]
            img_filename = filename.substr(0, filename.find_last_of(".")-4) + ".img.gz";
            img_output.open(img_filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }
        else if (extension == ".hdr.gz") {
            hdr_output.open(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

            // Change [.hdr.gz] to [.img.gz]
            img_filename = filename.substr(0, filename.find_last_of(".")-7) + ".img.gz";
            img_output.open(img_filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }
        else {
            hdr_output.open(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

            // Change [.hdr] to [.img]
            img_filename = filename.substr(0, filename.find_last_of(".")-4) + ".img";
            img_output.open(img_filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }

        // Check if NifTI-1 header output stream file is open.
        if (!hdr_output.is_open()) {
            throw std::runtime_error(Logger::Error("Could not open header file to save voxelized image in "
                                                   "NifTI-1 image format [.hdr | .hdr.gz]. "
                                                   "Check output filename: " + filename).c_str());
        }

        // Check if NifTI-1 header output stream file is open.
        if (!img_output.is_open()) {
            throw std::runtime_error(Logger::Error("Could not open image file to save voxelized image in "
                                                   "NifTI-1 image format [.img | .img.gz]. "
                                                   "Check output filename: " + img_filename).c_str());
        }

        // Write uncompressed header binary data in stringstream.
        std::stringstream uncomp_hdr_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);
        uncomp_hdr_data.write(reinterpret_cast<char *>(&ni1_header), sizeof(ni1_header));

        // Write uncompressed image binary data in stringstream.
        std::stringstream uncomp_img_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);
        uncomp_img_data.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

        // Write in files compressed mode.
        if (compression == ImageCompression::ON || extension == ".hdr.gz") {

            // Sting streams to store the compressed header and image binary data.
            std::stringstream comp_hdr_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);
            std::stringstream comp_img_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

            // Compressed file writing scope.
            {
                // Create compressed output header buffer.
                boost::iostreams::filtering_streambuf<boost::iostreams::input> comp_hdr_buffer;
                comp_hdr_buffer.push(boost::iostreams::gzip_compressor());
                comp_hdr_buffer.push(uncomp_hdr_data);

                // Copy compressed output header buffer in compressed binary data string stream.
                boost::iostreams::copy(comp_hdr_buffer, comp_hdr_data);

                // Flush the data from the compressor.
                boost::iostreams::close(comp_hdr_buffer);

                // Write the compressed header binary data.
                hdr_output.write(comp_hdr_data.str().c_str(), comp_hdr_data.str().length());

                // Create compressed output image buffer.
                boost::iostreams::filtering_streambuf<boost::iostreams::input> comp_img_buffer;
                comp_img_buffer.push(boost::iostreams::gzip_compressor());
                comp_img_buffer.push(uncomp_img_data);

                // Copy compressed output image buffer in compressed binary data string stream.
                boost::iostreams::copy(comp_img_buffer, comp_img_data);

                // Flush the data from the compressor.
                boost::iostreams::close(comp_img_data);

                // Write the compressed image binary data.
                img_output.write(comp_img_data.str().c_str(), comp_img_data.str().length());

            } // End of Compressed file writing scope.

        }
        else { // Write files in uncompressed mode.
            hdr_output.write(uncomp_hdr_data.str().c_str(), uncomp_hdr_data.str().length());
            img_output.write(uncomp_img_data.str().c_str(), uncomp_img_data.str().length());
        }

        // Close the output files.
        hdr_output.close();
        img_output.close();

    } // Write output.

}


std::ostream & operator << (std::ostream &out, const Nifti1Header &header)
{
    out << "header size: " << header.hdr_size_ << std::endl;
    out << "data type: " << header.data_type_ << std::endl;
    out << "db_name: " << header.db_name_ << std::endl;
    out << "extents: " << header.extents_ << std::endl;
    out << "session_error: " << header.session_error_ << std::endl;
    out << "regular: " << header.regular_ << std::endl;
    out << "dim info: " << header.dim_info_ << std::endl;
    for (int i=0; i != 8; ++i)
        out << "DIM" << i << " " << header.dim_[i] << std::endl;
    out << "intent p1: " << header.intent_p1_ << std::endl;
    out << "intent p2: " << header.intent_p2_ << std::endl;
    out << "intent p3: " << header.intent_p3_ << std::endl;
    out << "intent code: " << header.intent_code_ << std::endl;
    out << "datatype: " << header.datatype_ << std::endl;
    out << "bitpix: " << header.bitpix_ << std::endl;
    out << "slice start: " << header.slice_start_ << std::endl;
    for (int i=0; i != 8; ++i)
        out << "PIX DIM[" << i << "]: " << header.pixdim_[i] << std::endl;
    out << "vox Origin: " << header.vox_offset_ << std::endl;
    out << "scale slope: " << header.scl_slope_ << std::endl;
    out << "scale inter: " << header.scl_inter_ << std::endl;
    out << "slice end: " << header.slice_end_ << std::endl;
    out << "slice code: " << header.slice_code_ << std::endl;
    out << "xyzt units: " << header.xyzt_units_ << std::endl;
    out << "cal max: " << header.cal_max_ << std::endl;
    out << "cal min: " << header.cal_min_ << std::endl;
    out << "slice duration: " << header.slice_duration_ << std::endl;
    out << "time Origin: " << header.toffset_ << std::endl;
    out << "glmax: " << header.glmax_ << std::endl;
    out << "glmin: " << header.glmin_ << std::endl;
    out << "descrption: " << header.descrip_ << std::endl;
    out << "aux file: " << header.aux_file_ << std::endl;
    out << "qform code: " << header.qform_code_ << std::endl;
    out << "sform code: " << header.sform_code_ << std::endl;
    out << "quaternion b: " << header.quatern_b_ << std::endl;
    out << "quaternion c: " << header.quatern_c_ << std::endl;
    out << "quaternion d: " << header.quatern_d_ << std::endl;
    out << "Origin x: " << header.qoffset_x_ << std::endl;
    out << "Origin y: " << header.qoffset_y_ << std::endl;
    out << "Origin z: " << header.qoffset_z_ << std::endl;
    for (int i=0; i != 4; ++i)
        out << "srow x[" << i << "]: " << header.srow_x_[i] << std::endl;
    for (int i=0; i != 4; ++i)
        out << "srow y[" << i << "]: " << header.srow_y_[i] << std::endl;
    for (int i=0; i != 4; ++i)
        out << "srow z[" << i << "]: " << header.srow_z_[i] << std::endl;
    out << "intent name: " << header.intent_name_ << std::endl;
    out << "magic: " << header.magic_ << std::endl;

    return out;
}


} //end of namespace IMP

#endif //IMP_ENGINE_IMAGE_IO_NII_IO_TPP_
