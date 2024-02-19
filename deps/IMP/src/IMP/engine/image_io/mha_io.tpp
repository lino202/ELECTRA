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


#ifndef IMP_ENGINE_IMAGE_IO_MHA_IO_TPP_
#define IMP_ENGINE_IMAGE_IO_MHA_IO_TPP_


#include "IMP/engine/image_io/mha_io.hpp"


namespace IMP {


MhaIO::MhaIO()
{}


MhaIO::~MhaIO()
{}


template<typename VOXELTYPE>
void MhaIO::Load(const std::string &filename, VoxelImage<VOXELTYPE> &image)
{
    // Check file extension.
    std::string extension = filename.substr(filename.find_last_of("."));
    if (extension != ".mha") {
        std::string error = "[IMP ERROR] Cannot load MetaImage of unknown format. Expected format: [.mhd]";
        throw std::invalid_argument(error.c_str());
    }

    // Open MetaImage image file.
    std::ifstream mha_file(filename.c_str(), std::ios::in | std::ios::binary);

    if (!mha_file.is_open()) {
        std::string error = "[IMP ERROR] Cannot open the MetaImage file (.mha). Check given path.";
        throw std::runtime_error(error.c_str());
    }

    // Count the characters of the MetaImage header.
    std::string line = "";
    int chars_num = 0;
    while(std::getline(mha_file, line)) {

        // Transform header line to lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Count the length of the line plus the escaped '\n' delimiter.
        chars_num += line.length()+1;

        // Break the loop when the last line of the header is counted.
        if (line.find("elementdatafile") != std::string::npos) { break; }

    }

    // Return to the beginning of the image file.
    mha_file.clear();
    mha_file.seekg(0, mha_file.beg);

    // Extract MetaImage header information.
    std::string header_lines(chars_num, '\0');
    mha_file.read(&header_lines[0], chars_num);

    // Convert header to string stream for processing.
    std::stringstream header_stream(header_lines);

    //Header attributes.
    std::string object_type = "", num_dimensions = "", binary_data = "";
    std::string big_endian_order = "";
    std::string compressed_data = "";
    IMP::Vec<3, int> dimensions({-1, -1, -1});
    IMP::Vec<3, double> Origin({0., 0., 0.}), spacing({-1., -1., -1.});
    std::string element_type = "", element_data_file = "";

    // Read the header of the MetaImage.
    line = "";
    while(std::getline(header_stream, line)) {
        // Transform header line in lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Get equal sign position.
        auto equal_pos = line.find("=");

        // Replace equal sign if found with blank space.
        if (equal_pos != std::string::npos){ line.replace(equal_pos, 1," "); }

        // Get image attributes from the header file.
        if (line.find("objecttype") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, object_type); }

        else if (line.find("ndims") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, num_dimensions); }

        else if (line.find("binarydata") != std::string::npos && line.find("byteordermsb") == std::string::npos) {
            this->ExtractAttributeFromHeader<std::string>(line, binary_data);
        }

        else if (line.find("binarydatabyteordermsb") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, big_endian_order); }

        else if (line.find("compresseddata") != std::string::npos && line.find("size") == std::string::npos) {
            this->ExtractAttributeFromHeader<std::string>(line, compressed_data);
        }

        else if (line.find("Origin") != std::string::npos) { this->ExtractVecAttributeFromHeader<double>(line, Origin); }

        else if (line.find("elementspacing") != std::string::npos) { this->ExtractVecAttributeFromHeader<double>(line, spacing); }

        else if (line.find("dimsize") != std::string::npos) { this->ExtractVecAttributeFromHeader<int>(line, dimensions);}

        else if (line.find("elementtype") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, element_type); }

        else if (line.find("elementdatafile") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, element_data_file); }

        else {}

    }

    //Check for loading errors.
    if (object_type != "image") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mha). Expected object type is image.";
        throw std::runtime_error(error.c_str());
    }
    if (num_dimensions != "3") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mha) with more than 3 Dimensions.";
        throw std::runtime_error(error.c_str());
    }
    if (binary_data != "true") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mha). Expected binary image type.";
        throw std::runtime_error(error.c_str());
    }
    if (big_endian_order != "false") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mha). Expected low endianness.";
        throw std::runtime_error(error.c_str());
    }
    if (spacing[0] == -1 || spacing[1] == -1. || spacing[2] == -1.) {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mha). Spacing was not initialized correctly.";
        throw std::runtime_error(error.c_str());
    }
    if (dimensions[0] == -1 || dimensions[1] == -1 || dimensions[2] == -1) {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mha). Image dimensions were not initialized correctly.";
        throw std::runtime_error(error.c_str());
    }

    // Initialize image.
    image.SetDimensions(dimensions[0], dimensions[1], dimensions[2]);
    image.SetSpacing(spacing[0], spacing[1], spacing[2]);
    image.SetOrigin(Origin[0], Origin[1], Origin[2]);
    image.SetNumVoxels(dimensions[0]*dimensions[1]*dimensions[2]);
    image.SetMemorySize(static_cast<int>(sizeof(VOXELTYPE))*image.NumVoxels());

    // Set image data compression.
    if (compressed_data == "true") { image.SetCompression(ImageCompression::ON); }
    else { image.SetCompression(ImageCompression::OFF); }

    // Initialize the data type of the image.
    if (element_type == "met_char") { image.SetDataType(ImageDataType::char_type); }
    else if (element_type == "met_uchar") { image.SetDataType(ImageDataType::uchar_type); }
    else if (element_type == "met_short") { image.SetDataType(ImageDataType::short_int_type); }
    else if (element_type == "met_ushort") { image.SetDataType(ImageDataType::ushort_int_type); }
    else if (element_type == "met_int") { image.SetDataType(ImageDataType::int_type); }
    else if (element_type == "met_uint") { image.SetDataType(ImageDataType::uint_type); }
    else if (element_type == "met_long") { image.SetDataType(ImageDataType::long_int_type); }
    else if (element_type == "met_ulong") { image.SetDataType(ImageDataType::ulong_int_type); }
    else if (element_type == "met_float") { image.SetDataType(ImageDataType::float_type); }
    else if (element_type == "met_double") { image.SetDataType(ImageDataType::double_type); }
    else if (element_type == "met_bool") { image.SetDataType(ImageDataType::bool_type); }
    else { std::string error = "[IMP ERROR] Can not load MetaImage (.mha) of unknown element type.";
           throw std::runtime_error(error.c_str());
    }

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
         !(image.DataType() == ImageDataType::double_type && sizeof(VOXELTYPE) == sizeof(double)) &&
         !(image.DataType() == ImageDataType::bool_type && sizeof(VOXELTYPE) == sizeof(bool))
         ) {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mha) of incompatible data type.";
        throw std::runtime_error(error.c_str());
    }

    // Reset image data container.
    image.EditData().clear();
    image.EditData().assign(image.NumVoxels(), 0);

    if (compressed_data == "true") {
        // Compressed file reading scope.
        {
            // Read compressed binary data.
            boost::iostreams::filtering_istream input;
            input.push(boost::iostreams::zlib_decompressor());
            input.push(mha_file);

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
    else {
        // Create binary data buffer.
        VOXELTYPE *buffer = new VOXELTYPE [image.NumVoxels()];

        // Read Metaimage binary data.
        mha_file.read((char *)buffer, image.MemorySize());

        // Store the binary data in the voxelized image's data vector.
        for(auto i = 0; i != image.NumVoxels(); ++i) {
            image.EditData()[i] = static_cast<VOXELTYPE>(buffer[i]);
        }

        // Clean-up buffer allocated memory.
        delete[] buffer;
        buffer = nullptr;

    }


    // Close Metaimage file.
    mha_file.close();

}


template<typename VOXELTYPE>
void MhaIO::Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression)
{
    // MetaImage data type.
    std::string elem_type = "";

    switch ( image.DataType() ){
        case ImageDataType::char_type :
            elem_type = "MET_CHAR";
            break;
        case ImageDataType::uchar_type :
            elem_type = "MET_UCHAR";
            break;
        case ImageDataType::short_int_type :
            elem_type = "MET_SHORT";
            break;
        case ImageDataType::ushort_int_type :
            elem_type = "MET_USHORT";
            break;
        case ImageDataType::int_type :
            elem_type = "MET_INT";
            break;
        case ImageDataType::uint_type :
            elem_type = "MET_UINT";
            break;
        case ImageDataType::long_int_type :
            elem_type = "MET_LONG";
            break;
        case ImageDataType::ulong_int_type :
            elem_type = "MET_ULONG";
            break;
        case ImageDataType::float_type :
            elem_type = "MET_FLOAT";
            break;
        case ImageDataType::double_type :
            elem_type = "MET_DOUBLE";
            break;
        case ImageDataType::bool_type :
            elem_type = "MET_BOOL";
            break;
        default:
            std::string error = "[IMP ERROR] Cannot save image of unknown data type.";
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

    // Create the header file of the MetaImage according to the voxelized image attributes.
    std::string header = "ObjectType = Image\n";
    header += "NDims = 3\n";
    header += "BinaryData = True\n";
    header += "BinaryDataByteOrderMSB = False\n";

    // Get the position in the header where compression info will be inserted later.
    std::size_t comp_info_pos = header.length();

    header += "TransformMatrix = 1 0 0 0 1 0 0 0 1\n";
    header += "Origin = " + std::to_string(image.Origin()[0]) + " ";
    header += std::to_string(image.Origin()[1]) + " ";
    header += std::to_string(image.Origin()[2]) + "\n";
    header += "CenterOfRotation = 0 0 0\n";
    header += "AnatomicalOrientation = RAI\n";
    header += "ElementSpacing = " + std::to_string(image.Spacing()[0]) + " ";
    header += std::to_string(image.Spacing()[1]) + " ";
    header += std::to_string(image.Spacing()[2]) + "\n";
    header += "DimSize = " + std::to_string(image.Dimensions()[0]) + " ";
    header += std::to_string(image.Dimensions()[1]) + " ";
    header += std::to_string(image.Dimensions()[2]) + "\n";
    header += "ElementType = " + elem_type + "\n";
    header += "ElementDataFile = LOCAL\n";

    // Open the MetaImage output stream file.
    std::ofstream mha_output(filename, std::ios::out | std::ios::binary | std::ios::trunc);

    if (!mha_output.is_open()) {
        std::string error = "[IMP ERROR] Could not open file to save voxelized image in MetaImage format (.mha).";
        throw std::invalid_argument(error.c_str());
    }

    // Write the image binary data.
    if (compression == ImageCompression::ON) {

        // Convert image binary data to sting stream.
        std::stringstream uncompressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);
        uncompressed_data.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

        // Sting stream to store the compressed binary data.
        std::stringstream compressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

        // Set the image's compressed state ON.
        image.SetCompression(ImageCompression::ON);

        // Compressed file writing scope.
        {
            // Create compressed output buffer.
            boost::iostreams::filtering_streambuf<boost::iostreams::input> comp_buffer;
            comp_buffer.push(boost::iostreams::zlib_compressor());
            comp_buffer.push(uncompressed_data);

            // Copy compressed output buffer in compressed binary data string stream.
            boost::iostreams::copy(comp_buffer, compressed_data);

            // Flush the data from the compressor.
            boost::iostreams::close(comp_buffer);

            // Set the number of the compressed bytes of the image.
            image.SetCompressedBytes(compressed_data.str().length());

        } // end of Compressed file writing scope.

        // The compression info.
        std::string comp_info = "CompressedData = True\n";
        comp_info += "CompressedDataSize = " + std::to_string(image.CompressedBytes()) + "\n";

        // Add compression info in the header.
        header.insert(comp_info_pos, comp_info);

        // Write the header data.
        mha_output.write(header.c_str(), header.length());

        // Write the compressed image binary data.
        mha_output.write(compressed_data.str().c_str(), compressed_data.str().length());

    }
    else {
        // The compression info.
        std::string comp_info = "CompressedData = False\n";

        // Add compression info in the header.
        header.insert(comp_info_pos, comp_info);

        // Write the header data.
        mha_output.write(header.c_str(), header.length());

        // Write the image binary data.
        mha_output.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

    }

    // Close MetaImage output stream file.
    mha_output.close();

}


template <typename DATATYPE>
void MhaIO::ExtractAttributeFromHeader(const std::string &line, DATATYPE &attribute)
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


template <typename DATATYPE>
void MhaIO::ExtractVecAttributeFromHeader(const std::string &line, Vec<3, DATATYPE> &vec)
{
    // Convert line string to stringstream.
    std::stringstream ss(line);

    // First part of header's line.
    std::string token = "";

    // Attribute values to be extracted.
    DATATYPE x = 0; DATATYPE y = 0; DATATYPE z = 0;

    // Get the attribute values by processing the stringstream.
    if (!(ss >> token >> x >> y >> z)) {
        std::string error = "[IMP ERROR] Cannot extract 3-dimensional attribute from the MetaImage header (.mha) from line: " + line;
        throw std::runtime_error(error.c_str());
    }

    // Populate the vector with the attributes.
    vec.Set({x, y, z});

}


} //end of namespace IMP


#endif //IMP_ENGINE_IMAGE_IO_MHA_IO_TPP_
