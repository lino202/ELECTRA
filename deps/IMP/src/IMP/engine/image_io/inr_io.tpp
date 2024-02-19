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


#ifndef IMP_ENGINE_IMAGE_IO_INR_IO_TPP_
#define IMP_ENGINE_IMAGE_IO_INR_IO_TPP_


#include "IMP/engine/image_io/inr_io.hpp"


namespace IMP {


InrIO::InrIO()
{}


InrIO::~InrIO()
{}


template<typename VOXELTYPE>
void InrIO::Load(const std::string &filename, VoxelImage<VOXELTYPE> &image)
{
    // The Inria image string stream.
    std::stringstream inr_sstream(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

    // Open Inria image according to the extension.
    std::string extension = filename.substr(filename.find_last_of("."));
    if (extension == ".inr") {
        // Open Inria image file.
        std::ifstream inr_file(filename.c_str(), std::ios_base::in | std::ios_base::binary);

        // Check if file is open.
        if (!inr_file.is_open()) {
            std::string error = "[IMP ERROR] Cannot open the Inria image file (.inr). Check given path.";
            throw std::runtime_error(error.c_str());
        }

        // Read the image data.
        inr_sstream << inr_file.rdbuf();
        inr_file.close();
    }
    else {
        // Search for complete extension.
        extension = filename.substr(filename.find_last_of(".")-4);
        if (extension == ".inr.gz") {
            // Open Inria compressed image file.
            std::ifstream inr_comp_file(filename.c_str(), std::ios_base::in | std::ios_base::binary);

            // Check if file is open.
            if (!inr_comp_file.is_open()) {
                std::string error = "[IMP ERROR] Cannot open the Inria image file (.inr.gz). Check given path.";
                throw std::runtime_error(error.c_str());
            }

            // Image decompressing scope.
            {
                // Read compressed binary data.
                boost::iostreams::filtering_istream input;
                input.push(boost::iostreams::gzip_decompressor());
                input.push(inr_comp_file);
                boost::iostreams::copy(input, inr_sstream);
            } // end of Image decompressing scope.

            inr_comp_file.close();
        }
        else {
            std::string error = "[IMP ERROR] Cannot load Inria image of unknown format. Expected format: [.inr | .inr.gz]";
            throw std::invalid_argument(error.c_str());
        }
    }

    // Extract Inria image header information.
    std::string header_lines(256, '\0');
    inr_sstream.read(&header_lines[0], 256);

    // Convert header information to string stream.
    std::stringstream header_sstream(header_lines);
    //Header attributes.
    int x_dim = -1, y_dim = -1, z_dim = -1;
    int vector_dim = 0;
    std::string data_type = "";
    int voxel_bits = 0;
    std::string cpu = "";
    double x_spacing = -1., y_spacing = -1., z_spacing = -1.;
    double x_Origin = 0., y_Origin = 0., z_Origin = 0.;

    // Read the header of the Inria image.
    std::string line = "";
    while(std::getline(header_sstream, line)) {

        // Transform header stream line in lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Get equal sign position.
        auto equal_pos = line.find("=");

        // Replace equal sign if found with blank space.
        if (equal_pos != std::string::npos){ line.replace(equal_pos, 1," "); }

        // Get image attributes from the header file.
        if (line.find("xdim") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, x_dim); }

        if (line.find("ydim") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, y_dim); }

        if (line.find("zdim") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, z_dim); }

        if (line.find("vdim") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, vector_dim); }

        if (line.find("type") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, data_type); }

        if (line.find("pixsize") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, voxel_bits); }

        if (line.find("cpu") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, cpu); }

        if (line.find("vx") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, x_spacing); }

        if (line.find("vy") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, y_spacing); }

        if (line.find("vz") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, z_spacing); }

        if (line.find("tx") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, x_Origin); }

        if (line.find("ty") != std::string::npos && line.find("pe") == std::string::npos) { this->ExtractAttributeFromHeader<double>(line, y_Origin); }

        if (line.find("tz") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, z_Origin); }

    }

    // Check for loading errors.
    if (x_dim == -1 || y_dim == -1 || z_dim == -1) {
        std::string error = "[IMP ERROR] Cannot load Inria image (.inr). Image dimensions were not initialized correctly.";
        throw std::runtime_error(error.c_str());
    }
    if (vector_dim != 1) {
        std::string error = "[IMP ERROR] Cannot load Inria image (.inr). Only scalar images are currently supported.";
        throw std::runtime_error(error.c_str());
    }
    if ((voxel_bits/8) != sizeof(VOXELTYPE)) {
        std::string error = "[IMP ERROR] Cannot load Inria image (.inr) of incompatible data type.";
        throw std::runtime_error(error.c_str());
    }
    if (cpu != "decm") {
        std::string error = "[IMP ERROR] Cannot load Inria image (.inr). Expected low endianness.";
        throw std::runtime_error(error.c_str());
    }
    if (x_spacing == -1. || y_spacing == -1. || z_spacing == -1.) {
        std::string error = "[IMP ERROR] Cannot load Inria image (.inr). Spacing was not initialized correctly.";
        throw std::runtime_error(error.c_str());
    }

    // Initialize image.
    image.SetDimensions(x_dim, y_dim, z_dim);
    image.SetSpacing(x_spacing, y_spacing, z_spacing);
    image.SetOrigin(x_Origin, y_Origin, z_Origin);
    image.SetNumVoxels(x_dim*y_dim*z_dim);
    image.SetMemorySize(static_cast<int>(sizeof(VOXELTYPE))*image.NumVoxels());

    // Set image data compression.
    if (extension == ".inr.gz") { image.SetCompression(ImageCompression::ON); }
    else { image.SetCompression(ImageCompression::OFF); }

    // Initialize the data type of the image.
    if (data_type == "signed" && voxel_bits == 8) { image.SetDataType(ImageDataType::char_type); }
    else if (data_type == "unsigned" && voxel_bits == 8) { image.SetDataType(ImageDataType::uchar_type); }
    else if (data_type == "signed" && voxel_bits == 16) { image.SetDataType(ImageDataType::short_int_type); }
    else if (data_type == "unsigned" && voxel_bits == 16) { image.SetDataType(ImageDataType::ushort_int_type); }
    else if (data_type == "signed" && voxel_bits == 32) { image.SetDataType(ImageDataType::int_type); }
    else if (data_type == "unsigned" && voxel_bits == 32) { image.SetDataType(ImageDataType::uint_type); }
    else if (data_type == "signed" && voxel_bits == 64) { image.SetDataType(ImageDataType::long_int_type); }
    else if (data_type == "unsigned" && voxel_bits == 64) { image.SetDataType(ImageDataType::ulong_int_type); }
    else if (data_type == "float" && voxel_bits == 32) { image.SetDataType(ImageDataType::float_type); }
    else if (data_type == "float" && voxel_bits == 64) { image.SetDataType(ImageDataType::double_type); }
    else if (data_type == "bool" && voxel_bits == 8) { image.SetDataType(ImageDataType::bool_type); }
    else { std::string error = "[IMP ERROR] Cannot load Inria image (.inr) of unknown data type.";
           throw std::runtime_error(error.c_str());
    }

    // Reset image data container.
    image.EditData().clear();
    image.EditData().assign(image.NumVoxels(), 0);

    // Create binary data buffer.
    VOXELTYPE *buffer = new VOXELTYPE [image.NumVoxels()];

    // Read Inria image data.
    inr_sstream.read((char *)buffer, image.MemorySize());

    // Store the binary data in the voxelized image's data vector.
    for(auto i = 0; i != image.NumVoxels(); ++i) {
        image.EditData()[i] = static_cast<VOXELTYPE>(buffer[i]);
    }

    // Clean-up buffer allocated memory.
    delete[] buffer;
    buffer = nullptr;

}


template<typename VOXELTYPE>
void InrIO::Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression)
{
    // Inria image data type.
    std::string data_type = "";

    // Inria image number of bits per voxel.
    int voxel_bits = 0;

    // Set Inria image data type and voxel bits according to the data type of the given voxelized image.
    switch (image.DataType()) {
        case IMP::ImageDataType::char_type :
            data_type = "signed fixed";
            voxel_bits = 8;
            break;
        case IMP::ImageDataType::uchar_type :
            data_type = "unsigned fixed";
            voxel_bits = 8;
            break;
        case IMP::ImageDataType::short_int_type :
            data_type = "signed fixed";
            voxel_bits = 16;
            break;
        case IMP::ImageDataType::ushort_int_type :
            data_type = "unsigned fixed";
            voxel_bits = 16;
            break;
        case IMP::ImageDataType::int_type :
            data_type = "signed fixed";
            voxel_bits = 32;
            break;
        case IMP::ImageDataType::uint_type :
            data_type = "unsigned fixed";
            voxel_bits = 32;
            break;
        case IMP::ImageDataType::long_int_type :
            data_type = "signed fixed";
            voxel_bits = 64;
            break;
        case IMP::ImageDataType::ulong_int_type :
            data_type = "unsigned fixed";
            voxel_bits = 64;
            break;
        case IMP::ImageDataType::float_type :
            data_type = "float";
            voxel_bits = 32;
            break;
        case IMP::ImageDataType::double_type :
            data_type = "float";
            voxel_bits = 64;
            break;
        case IMP::ImageDataType::bool_type :
            data_type = "bool";
            voxel_bits = 8;
            break;
        default:
            std::string error = "[IMP ERROR] Cannot save image of unknown data type in Inria image format.";
            throw std::invalid_argument(error.c_str());
    }

    // Initialize the path of the output file.
    std::string path = "";

    // Position of the last slash in the output file's name.
    std::size_t last_slash = filename.find_last_of("/\\");

    // Get the path directory of the output file.
    if (last_slash != std::string::npos) {
        path = filename.substr(0, last_slash);
    }

    // Create the path's directory if it doesn't exist.
    boost::filesystem::path dir(path);
    if (!path.empty() && !boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
    }

    // Create the header file of the Inria image according to the voxelized image attributes.
    std::string header = "#INRIMAGE-4#{\n";
    header += "XDIM=" + std::to_string(image.Dimensions()[0]) + "\n";
    header += "YDIM=" + std::to_string(image.Dimensions()[1]) + "\n";
    header += "ZDIM=" + std::to_string(image.Dimensions()[2]) + "\n";
    header += "VDIM=1\n";
    header += "TYPE=" + data_type + "\n";
    header += "PIXSIZE=" + std::to_string(voxel_bits) + " bits\n";
    header += "CPU=decm\n";
    header += "VX=" + std::to_string(image.Spacing()[0]) + "\n";
    header += "VY=" + std::to_string(image.Spacing()[1]) + "\n";
    header += "VZ=" + std::to_string(image.Spacing()[2]) + "\n";
    header += "TX=" + std::to_string(image.Origin()[0]) + "\n";
    header += "TY=" + std::to_string(image.Origin()[1]) + "\n";
    header += "TZ=" + std::to_string(image.Origin()[2]) + "\n";

    // Fill with new line characters until header is 252 bytes long.
    header.append(252 - header.length(), '\n');

    // Add end of header. Header size now 256 bytes long.
    header += "##}\n";

    // Get Inria output file extension.
    std::string extension = filename.substr(filename.find_last_of("."));
    if (extension != ".inr") {
        // Search for complete extension.
        extension = filename.substr(filename.find_last_of(".")-4);
        if (extension != ".inr.gz") {
            std::string error = "[IMP ERROR] Cannot save Inria image to unsupported format. Supported format: [.inr | .inr.gz]";
            throw std::invalid_argument(error.c_str());
        }
    }

    // Open the Inria output stream file.
    std::ofstream inr_output;
    if (compression == ImageCompression::ON || extension == ".inr.gz") {

        // Open output stream file with complete compressed extension.
        if (extension == ".inr") {
            inr_output.open(filename + ".gz", std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }
        else {
            inr_output.open(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }

        // Check if Inria output stream file is open.
        if (!inr_output.is_open()) {
            std::string error = "[IMP ERROR] Could not open file to save voxelized image in "
                                "compressed Inria image format [.inr.gz].";
            throw std::runtime_error(error.c_str());
        }

        // Convert image binary data to sting stream.
        std::stringstream uncompressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);
        uncompressed_data.write(header.c_str(), header.length());
        uncompressed_data.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

        // Sting stream to store the compressed binary data.
        std::stringstream compressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

        // Set the image's compressed state ON.
        image.SetCompression(ImageCompression::ON);

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
            inr_output.write(compressed_data.str().c_str(), compressed_data.str().length());

            // Close Inria output stream file.
            inr_output.close();

        } // end of Compressed file writing scope.

    }
    else {
        inr_output.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);

        if (!inr_output.is_open()) {
            std::string error = "[IMP ERROR] Could not open file to save voxelized image in Inria image format (.inr).";
            throw std::runtime_error(error.c_str());
        }

         // Write the header data.
         inr_output.write(header.c_str(), header.length());

        // Write the image data.
        inr_output.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

        // Close Inria output stream file.
        inr_output.close();
    }

}


template <typename DATATYPE>
void InrIO::ExtractAttributeFromHeader(const std::string &line, DATATYPE &attribute)
{
    // Convert line string to stringstream.
    std::stringstream ss(line);

    // First part of header's line.
    std::string token = "";

    // Get the attribute by processing the stringstream.
    if (!(ss >> token >> attribute)) {
        std::string error = "[IMP ERROR] Can not extract attribute from the Inria image header (.inr) from line: " + line;
        throw std::runtime_error(error.c_str());
    }

}


} //end of namespace IMP


#endif //IMP_ENGINE_IMAGE_IO_INR_IO_TPP_
