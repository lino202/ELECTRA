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


#ifndef IMP_ENGINE_IMAGE_IO_MHD_IO_TPP_
#define IMP_ENGINE_IMAGE_IO_MHD_IO_TPP_


#include "IMP/engine/image_io/mhd_io.hpp"


namespace IMP {


MhdIO::MhdIO()
{}


MhdIO::~MhdIO()
{}


template<typename VOXELTYPE>
void MhdIO::Load(const std::string &filename, IMP::VoxelImage<VOXELTYPE> &image)
{
    //Check file extension.
    std::string extension = filename.substr(filename.find_last_of("."));
    if (extension != ".mhd") {
        std::string error = "[IMP ERROR] Cannot load metaimage of unknown format. Expected format: [.mhd]";
        throw std::invalid_argument(error.c_str());
    }

    //Open file.
    std::ifstream mhd_file(filename.c_str(), std::ios::in);

    if (!mhd_file.is_open()) {
        std::string error = "[IMP ERROR] Cannot open the metaimage file. Check given path.";
        throw std::runtime_error(error.c_str());
    }

    //Line of the file.
    std::string line = "";

    //Header attributes.
    std::string object_type = "", num_dimensions = "", binary_data = "";
    std::string big_endian_order = "", element_type = "";
    std::string compressed_data = "", element_data_file = "";
    IMP::Vec<3, int> dimensions({-1, -1, -1});
    IMP::Vec<3, double> Origin({0., 0., 0.}), spacing({-1., -1., -1.});

    //Read header.
    while (std::getline(mhd_file, line)) {

        // Transform .mhd file line in lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Get equal sign position.
        auto equal_pos = line.find("=");

        // Replace equal sign if found with blank space.
        if (equal_pos != std::string::npos){ line.replace(equal_pos, 1," "); }

        // Get image attributes from the header file.
        if (line.find("objecttype") != std::string::npos) { this->ExtractStringFromHeader(line, object_type); }

        if (line.find("ndims") != std::string::npos) { this->ExtractStringFromHeader(line, num_dimensions); }

        if (line.find("binarydata") != std::string::npos && line.find("byteordermsb") == std::string::npos) {
            this->ExtractStringFromHeader(line, binary_data);
        }

        if (line.find("binarydatabyteordermsb") != std::string::npos) { this->ExtractStringFromHeader(line, big_endian_order); }

        if (line.find("compresseddata") != std::string::npos && line.find("size") == std::string::npos) {
            this->ExtractStringFromHeader(line, compressed_data);
        }

        if (line.find("Origin") != std::string::npos) { this->ExtractVecFromHeader<double>(line, Origin); }

        if (line.find("elementspacing") != std::string::npos) { this->ExtractVecFromHeader<double>(line, spacing); }

        if (line.find("dimsize") != std::string::npos) { this->ExtractVecFromHeader<int>(line, dimensions);}

        if (line.find("elementtype") != std::string::npos) { this->ExtractStringFromHeader(line, element_type); }

        if (line.find("elementdatafile") != std::string::npos) { this->ExtractStringFromHeader(line, element_data_file); }
    }

    // Close .mhd file.
    mhd_file.close();

    //Check for loading errors.
    if (object_type != "image") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mhd). Expected object type is image.";
        throw std::runtime_error(error.c_str());
    }
    if (num_dimensions != "3") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mhd) with more than 3 Dimensions.";
        throw std::runtime_error(error.c_str());
    }
    if (binary_data != "true") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mhd). Expected binary image type.";
        throw std::runtime_error(error.c_str());
    }
    if (big_endian_order != "false") {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mhd). Expected low endianness.";
        throw std::runtime_error(error.c_str());
    }
    if (spacing[0] == -1 || spacing[1] == -1. || spacing[2] == -1.) {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mhd). Spacing was not initialized correctly.";
        throw std::runtime_error(error.c_str());
    }
    if (dimensions[0] == -1 || dimensions[1] == -1 || dimensions[2] == -1) {
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mhd). Image dimensions were not initialized correctly.";
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

    // Initialize the data type of the image. NOTE that element_type string is in lower case.
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
    else { std::string error = "[IMP ERROR] Can not load MetaImage (.mhd) of unknown element type.";
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
        std::string error = "[IMP ERROR] Cannot load MetaImage (.mhd) of incompatible data type.";
        throw std::runtime_error(error.c_str());
    }

    //Prepare raw file path.
    std::string ElementDataFileWithPath = filename;
    auto lastindex = ElementDataFileWithPath.find_last_of(".");
    ElementDataFileWithPath = ElementDataFileWithPath.substr(0, lastindex);

    // Load image's binary data.
    if (compressed_data == "true") {
        // Load compressed binary data.
        ElementDataFileWithPath += ".zraw";
        this->ReadZrawData(ElementDataFileWithPath.c_str(), image);
    }
    else {
        // Load normal binary data.
        ElementDataFileWithPath += ".raw";
        this->ReadRawData(ElementDataFileWithPath.c_str(), image);
    }

}


template<typename VOXELTYPE>
void MhdIO::Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression)
{
    // Voxelized image data type.
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
            std::string error = "[IMP ERROR] Cannot save MetaImage [.mhd] of unknown data type.";
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

    // Create the binary data output filepath.
    std::string  raw_filename = filename.substr(0, filename.find_last_of(".") );

    // Write the binary data in the output file.
    if (compression == ImageCompression::ON) {
        raw_filename += ".zraw";
        this->WriteZRawData(image, raw_filename);
    }
    else {
        raw_filename += ".raw";
        this->WriteRawData(image, raw_filename);
    }

    //Create the header file.
    std::ofstream output_mhd;
    output_mhd.open(filename.c_str());

    if (!output_mhd.is_open()) {
        std::string error = "[IMP ERROR] Could not open header file to save voxelized image in metaimage format (.mhd).";
        throw std::invalid_argument(error.c_str());
    }

    output_mhd << "ObjectType = Image\n";
    output_mhd << "NDims = 3\n";
    output_mhd << "BinaryData = True\n";
    output_mhd << "BinaryDataByteOrderMSB = False\n";

    if (compression == ImageCompression::ON) {
        output_mhd << "CompressedData = True\n";
        output_mhd << "CompressedDataSize = " << image.CompressedBytes() << std::endl;
    }
    else {
        output_mhd << "CompressedData = False\n";
    }

    output_mhd << "TransformMatrix = 1 0 0 0 1 0 0 0 1\n";

    output_mhd << "Origin = " << image.Origin()[0] << " "
               << image.Origin()[1] << " "
               << image.Origin()[2] << std::endl;

    output_mhd << "CenterOfRotation = 0 0 0\n";
    output_mhd << "AnatomicalOrientation = RAI\n";

    output_mhd << "ElementSpacing = " << image.Spacing()[0] << " "
               << image.Spacing()[1] << " "
               << image.Spacing()[2] << std::endl;

    output_mhd << "DimSize = " << image.Dimensions()[0] << " "
               << image.Dimensions()[1] << " "
               << image.Dimensions()[2] << std::endl;

    output_mhd << "ElementType = " << elem_type << std::endl;

    if (last_slash != std::string::npos) {
        std::string short_raw_filename = raw_filename.substr(last_slash);
        output_mhd << "ElementDataFile = ." << short_raw_filename << std::endl;
    }
    else { output_mhd << "ElementDataFile = ./" << raw_filename << std::endl; }

    output_mhd.close();

}


template<typename VOXELTYPE>
void MhdIO::ReadRawData(const std::string &raw_filename, IMP::VoxelImage<VOXELTYPE> &image)
{
    // Load binary file.
    std::ifstream binary_file(raw_filename.c_str(), std::ios::in | std::ios::binary);

    // Check if binary file is opened.
    if (!binary_file.is_open()) {
        std::string error = "[IMP ERROR] Cannot open binary metaimage file (.raw).";
        throw std::invalid_argument(error.c_str());
    }

    //Check file extension.
    std::string extension = raw_filename.substr(raw_filename.find_last_of("."));
    if (extension != ".raw") {
        std::string error = "[IMP ERROR] Cannot read binary metaimage of unknown format. Expected format: (.raw)";
        throw std::invalid_argument(error.c_str());
    }

    // Reset voxelized image's data to zero.
    image.EditData().clear();
    image.EditData().assign(image.NumVoxels(), 0);

    // Create binary data buffer.
    VOXELTYPE *buffer = new VOXELTYPE [image.NumVoxels()];

    // Read the binary data in the buffer.
    binary_file.read((char *)buffer, image.MemorySize());
    binary_file.close();

    // Store the binary data in the voxelized image's data vector.
    for(auto i = 0; i != image.NumVoxels(); ++i) {
        image.EditData()[i] = static_cast<VOXELTYPE>(buffer[i]);
    }

    // Clean-up buffer allocated memory.
    delete[] buffer;
    buffer = nullptr;
}


template<typename VOXELTYPE>
void MhdIO::ReadZrawData(const std::string &zraw_filename, VoxelImage<VOXELTYPE> &image)
{
    //Check file extension.
    std::string extension = zraw_filename.substr(zraw_filename.find_last_of("."));
    if (extension != ".zraw") {
        std::string error = "[IMP ERROR] Cannot read compressed binary metaimage of unknown compression. Expected compression: (.zraw)";
        throw std::invalid_argument(error.c_str());
    }

    // Check if image is initialized.
    if (!image.IsInitialized()) {
        std::string error = "[IMP ERROR] The voxelized image is not initialized. Set dimensions, number of voxels and memory size.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the image data.
    image.EditData().clear();
    image.EditData().reserve(image.NumVoxels());

    // Compressed file reading scope.
    {
        // Read compressed file.
        std::ifstream file_in(zraw_filename, std::ios::in | std::ios_base::binary);
        boost::iostreams::filtering_istream input;
        input.push(boost::iostreams::zlib_decompressor());
        input.push(file_in);

        // Get character pointer to the image data vector.
        char *data_ptr = reinterpret_cast<char*>(image.EditData().data());

        // Assign data to the image data vector byte by byte.
        for (int i = 0; i != image.MemorySize(); ++i) {
            input.read(&data_ptr[i], 1);
        }

        // Close compressed binary file.
        file_in.close();

        // Release character pointer.
        data_ptr = nullptr;

    } // End of compressed file reading scope.

}


template<typename VOXELTYPE>
void MhdIO::WriteRawData(VoxelImage<VOXELTYPE> &image, const std::string &raw_filename)
{
    //Create the binary file.
    std::ofstream output_raw(raw_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

    if (!output_raw.is_open()) {
        std::string error = "[IMP ERROR] Could not open binary file to save voxelized image in metaimage format (.mhd).";
        throw std::invalid_argument(error.c_str());
    }

    // Write binary data.
    output_raw.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

    // Close binary file.
    output_raw.close();

}


template<typename VOXELTYPE>
void MhdIO::WriteZRawData(VoxelImage<VOXELTYPE> &image, const std::string &zraw_filename)
{

    // Create the compressed binary file.
    std::ofstream output_zraw(zraw_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

    // Check if compressed binary file can be opened.
    if (!output_zraw.is_open()) {
        std::string error = "[IMP ERROR] Could not open binary file to save voxelized image in metaimage format (.mhd).";
        throw std::invalid_argument(error.c_str());
    }

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
        comp_buffer.push(boost::iostreams::counter());
        comp_buffer.push(boost::iostreams::zlib_compressor());
        comp_buffer.push(uncompressed_data);

        // Copy compressed output buffer in compressed binary file.
        boost::iostreams::copy(comp_buffer, compressed_data);

        // Flush the data from the compressor.
        boost::iostreams::close(comp_buffer);

        // Set the number of the compressed bytes of the image.
        image.SetCompressedBytes(compressed_data.str().length());

    } // end of Compressed file writing scope.

    // Write the compressed image binary data.
    output_zraw.write(compressed_data.str().c_str(), compressed_data.str().length());

    // Close compressed binary file.
    output_zraw.close();

}


void MhdIO::ExtractStringFromHeader(const std::string &line, std::string &attribute)
{    
    // Convert line string to stringstream.
    std::stringstream ss(line);

    // First part of header's line.
    std::string token = "";

    // Get the attribute by processing the stringstream.
    if (!(ss >> token >> attribute)) {
        std::string error = "[IMP ERROR] Can not extract string from the .mhd header.";
        throw std::runtime_error(error.c_str());
    }
}


template <typename DATATYPE>
void MhdIO::ExtractVecFromHeader(const std::string &line, Vec<3, DATATYPE> &vec)
{
    // Convert line string to stringstream.
    std::stringstream ss(line);

    // First part of header's line.
    std::string token = "";

    // Values to be extracted.
    DATATYPE x = 0; DATATYPE y = 0; DATATYPE z = 0;

    // Get the values by processing the stringstream.
    if (!(ss >> token >> x >> y >> z)) {
        std::string error = "[IMP ERROR] Can not extract x, y, z values from the .mhd header.";
        throw std::runtime_error(error.c_str());
    }

    // Populate the vector with the values.
    vec.Set({x, y, z});

}


} //end of namespace IMP


#endif //IMP_ENGINE_IMAGE_IO_MHD_IO_TPP_
