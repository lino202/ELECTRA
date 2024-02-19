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


#ifndef IMP_ENGINE_IMAGE_IO_H33_IO_TPP_
#define IMP_ENGINE_IMAGE_IO_H33_IO_TPP_


#include "IMP/engine/image_io/h33_io.hpp"


namespace IMP {


H33IO::H33IO()
{}


H33IO::~H33IO()
{}


template<typename VOXELTYPE>
void H33IO::Load(const std::string &filename, VoxelImage<VOXELTYPE> &image)
{
    if (filename.empty()) {
        std::string error = "[IMP ERROR] Empty filename given during loading of interfile image (.h33).";
        throw std::invalid_argument(error.c_str());
    }

    //Check file extension.
    std::string extension = filename.substr(filename.find_last_of("."));
    if (extension != ".h33") {
        std::string error = "[IMP ERROR] Cannot load interfile image of unknown format. Expected format: (.h33)";
        throw std::invalid_argument(error.c_str());
    }

    //Open file.
    std::ifstream h33_file(filename.c_str(), std::ios::in);

    if (!h33_file.is_open()) {
        std::string error = "[IMP ERROR] Could not open the file: " + filename + " Check given path.";
        throw std::runtime_error(error.c_str());
    }

    //Line of the file.
    std::string line = "";

    //Header attributes.
    int data_Origin_in_bytes = 0;
    std::string data_compression = "", data_encode = "", imagedata_byte_order = "";
    int x_dim = -1, y_dim = -1, z_dim = -1;
    std::string number_format = "";
    int number_of_bytes_per_pixel = 0;
    double x_spacing = -1., y_spacing = -1., z_spacing = -1.;

    while (std::getline(h33_file, line)) {

        // Transform .h33 file line in lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Process line containing attribute.
        if (line.find('=') != std::string::npos){
            // Get equal sign position.
            auto equal_pos = line.find('=');

            // Replace equal sign with blank space.
            line.replace(equal_pos, 1," ");

            // Remove all previous blank spaces before equal sign.
            line.erase(std::remove(line.begin(), line.begin()+equal_pos, ' '), line.begin()+equal_pos);

        }

        if (line.find("dataOrigininbytes") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, data_Origin_in_bytes); }

        if (line.find("datacompression") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, data_compression); }

        if (line.find("dataencode") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, data_encode); }

        if (line.find("imagedatabyteorder") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, imagedata_byte_order); }

        if (line.find("matrixsize[1]") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, x_dim); }

        if (line.find("matrixsize[2]") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, y_dim); }

        if (line.find("numberformat") != std::string::npos) { this->ExtractAttributeFromHeader<std::string>(line, number_format); }

        if (line.find("numberofbytesperpixel") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, number_of_bytes_per_pixel); }

        if (line.find("scalingfactor(mm/pixel)[1]") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, x_spacing); }

        if (line.find("scalingfactor(mm/pixel)[2]") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, y_spacing); }

        if (line.find("numberofslices") != std::string::npos) { this->ExtractAttributeFromHeader<int>(line, z_dim); }

        if (line.find("slicethickness(pixels)") != std::string::npos) { this->ExtractAttributeFromHeader<double>(line, z_spacing); }
    }

    // Close .h33 file.
    h33_file.close();

    //Check for loading errors.
    if (data_Origin_in_bytes != 0) {
        std::string error = "[IMP ERROR] Cannot load Interfile (.h33) voxelized image with non-zero data Origin.";
        throw std::invalid_argument(error.c_str());
    }
    if ((data_compression != "none") && (data_compression != "gzip")) {
        std::string error = "[IMP ERROR] Cannot load .h33 voxelized image of unknown compression. [Supported compression: none, gzip]";
        throw std::invalid_argument(error.c_str());
    }
    if (data_encode != "none") {
        std::string error = "[IMP ERROR] Cannot load .h33 voxelized image with encoded data. [Supported encoding: none]";
        throw std::invalid_argument(error.c_str());
    }
    if (imagedata_byte_order != "littleendian") {
        std::string error = "[IMP ERROR] Cannot load .h33 voxelized image of unknown endianness. [Supported endianness: littleendian]";
        throw std::invalid_argument(error.c_str());
    }
    if (x_dim == -1 || y_dim == -1 || z_dim == -1) {
        std::string error = "[IMP ERROR] Cannot load .h33 voxelized image of unknown dimensions.";
        throw std::invalid_argument(error.c_str());
    }
    if (x_spacing == -1. || y_spacing == -1. || z_spacing == -1.) {
        std::string error = "[IMP ERROR] Cannot load .h33 voxelized image of unknown spacing.";
        throw std::invalid_argument(error.c_str());
    }
    if (number_of_bytes_per_pixel != static_cast<int>(sizeof(VOXELTYPE)) ) {
        std::string error = "[IMP ERROR] Cannot load Interfile image (.h33) of incompatible data type.";
        throw std::invalid_argument(error.c_str());
    }

    // Initialize voxelized image.
    image.SetDimensions(x_dim, y_dim, z_dim);
    image.SetSpacing(x_spacing, y_spacing, z_spacing);
    image.SetOrigin(0., 0., 0.);
    image.SetNumVoxels(x_dim*y_dim*z_dim);
    image.SetMemorySize(static_cast<int>(sizeof(VOXELTYPE))*image.NumVoxels());

    // Set image data compression.
    if (data_compression == "gzip") { image.SetCompression(IMP::ImageCompression::ON); }
    else { image.SetCompression(IMP::ImageCompression::OFF); }

    // Initialize the data type of the image.
    if ((number_format == "signed") && (number_of_bytes_per_pixel == 1)) {
            image.SetDataType(IMP::ImageDataType::char_type);
    }
    else if ((number_format == "unsigned") && (number_of_bytes_per_pixel == 1)) {
        image.SetDataType(IMP::ImageDataType::uchar_type);
    }
    else if ((number_format == "signed") && (number_of_bytes_per_pixel == 2)) {
        image.SetDataType(IMP::ImageDataType::short_int_type);
    }
    else if ((number_format == "unsigned") && (number_of_bytes_per_pixel == 2)) {
        image.SetDataType(IMP::ImageDataType::ushort_int_type);
    }
    else if ((number_format == "signed") && (number_of_bytes_per_pixel == 4)) {
        image.SetDataType(IMP::ImageDataType::int_type);
    }
    else if ((number_format == "unsigned") && (number_of_bytes_per_pixel == 4)) {
        image.SetDataType(IMP::ImageDataType::uint_type);
    }
    else if ((number_format == "signed") && (number_of_bytes_per_pixel == 8)) {
        image.SetDataType(IMP::ImageDataType::long_int_type);
    }
    else if ((number_format == "unsigned") && (number_of_bytes_per_pixel == 8)) {
        image.SetDataType(IMP::ImageDataType::ulong_int_type);
    }
    else if ((number_format == "float") && (number_of_bytes_per_pixel == 4)) {
        image.SetDataType(IMP::ImageDataType::float_type);
    }
    else if (((number_format == "float") || (number_format == "double")) && (number_of_bytes_per_pixel == 8)) {
        image.SetDataType(IMP::ImageDataType::double_type);
    }
    else if ((number_format == "boolean") && (number_of_bytes_per_pixel == 1)) {
        image.SetDataType(IMP::ImageDataType::bool_type);
    }
    else {
        std::string error = "[IMP ERROR] Can not load .h33 voxelized image of unknown data type.";
        throw std::invalid_argument(error.c_str());
    }

    //Convert .h33 file path to binary .i33.
    std::string ElementDataFileWithPath = filename;
    auto lastindex = ElementDataFileWithPath.find_last_of(".");
    ElementDataFileWithPath = ElementDataFileWithPath.substr(0, lastindex);

    // Load image's binary data.
    if (data_compression == "gzip") {
        // Load compressed binary data.
        ElementDataFileWithPath += ".i33.gz";
        // this->ReadI33GzData(ElementDataFileWithPath.c_str(), image);
    }
    else {
        // Load normal binary data.
        ElementDataFileWithPath += ".i33";
        this->ReadI33Data(ElementDataFileWithPath.c_str(), image);
    }

}


template<typename VOXELTYPE>
void H33IO::Save(VoxelImage<VOXELTYPE> &image, const std::string &filename, const ImageCompression &compression)
{
    // Voxelized image data type and number of bytes per pixel.
    std::string number_format = "";
    int number_of_bytes_per_pixel = 0;

    switch ( image.DataType() ){
        case IMP::ImageDataType::char_type :
            number_format = "signed character";
            number_of_bytes_per_pixel = 1;
            break;
        case IMP::ImageDataType::uchar_type :
            number_format = "unsigned character";
            number_of_bytes_per_pixel = 1;
            break;
        case IMP::ImageDataType::short_int_type :
            number_format = "signed integer";
            number_of_bytes_per_pixel = 2;
            break;
        case IMP::ImageDataType::ushort_int_type :
            number_format = "unsigned integer";
            number_of_bytes_per_pixel = 2;
            break;
        case IMP::ImageDataType::int_type :
            number_format = "signed integer";
            number_of_bytes_per_pixel = 4;
            break;
        case IMP::ImageDataType::uint_type :
            number_format = "unsigned integer";
            number_of_bytes_per_pixel = 4;
            break;
        case IMP::ImageDataType::long_int_type :
            number_format = "signed integer";
            number_of_bytes_per_pixel = 8;
            break;
        case IMP::ImageDataType::ulong_int_type :
            number_format = "unsigned integer";
            number_of_bytes_per_pixel = 8;
            break;
        case IMP::ImageDataType::float_type :
            number_format = "float";
            number_of_bytes_per_pixel = 4;
            break;
        case IMP::ImageDataType::double_type :
            number_format = "float";
            number_of_bytes_per_pixel = 8;
            break;
        case IMP::ImageDataType::bool_type :
            number_format = "unsigned character";
            number_of_bytes_per_pixel = 1;
            break;
        default:
            std::string error = "[IMP ERROR] Cannot save image of unsupported data type in interfile format (.h33).";
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
    std::string  i33_filename = filename.substr(0, filename.find_last_of(".") );

    // Write the binary data in the output file.
    if (compression == ImageCompression::ON) {
        i33_filename += ".i33.gz";
        // this->WriteI33GzData(image, i33_filename);
    }
    else {
        i33_filename += ".i33";
        this->WriteI33Data(image, i33_filename);
    }

    //Create the header file.
    std::ofstream output_h33(filename.c_str(), std::ios::out | std::ios::trunc );

    if (!output_h33.is_open()) {
        std::string error = "[IMP ERROR] Could not open header file to save voxelized image in interfile format (.h33).";
        throw std::invalid_argument(error.c_str());
    }

    output_h33 << "!INTERFILE :=\n";
    output_h33 << "!imaging modality := xx\n";
    output_h33 << "!originating system := xx\n";
    output_h33 << "!version of keys := 0.0\n";
    output_h33 << "date of keys := 0000:00:00\n";
    output_h33 << "conversion program := IMP Library\n";
    output_h33 << "program author := Konstantinos A. Mountris\n";
    output_h33 << "program version := 0.0.1\n";
    output_h33 << "program date := 00:00:00\n;\n";

    output_h33 << "!GENERAL DATA :=\n";
    output_h33 << "!data Origin in bytes := 0\n";

    if (last_slash != std::string::npos) {
        std::string short_i33_filename = i33_filename.substr(last_slash);
        output_h33 << "!name of data file := ." << short_i33_filename << std::endl;
    }
    else { output_h33 << "!name of data file := ./" << i33_filename << std::endl; }

    output_h33 << "patient name := xx\n";
    output_h33 << "!patient ID := 0000\n";
    output_h33 << "patient dob := 0000:00:00\n";
    output_h33 << "patient sex := Unknown\n";
    output_h33 << "!exam type := Unknown\n";

    if (compression == ImageCompression::ON) {
        output_h33 << "data compression := Gzip\n";
        output_h33 << "compressed data size := " << image.CompressedBytes() << std::endl;
    }
    else {
        output_h33 << "data compression := none\n";
    }

    output_h33 << "data encode := none\n;\n";

    output_h33 << "!GENERAL IMAGE DATA :=\n";
    output_h33 << "!type of data := Tomographic\n";
    output_h33 << "!total number of images := " << image.Dimensions()[2] << std::endl;
    output_h33 << "study date := 0000:00:00\n";
    output_h33 << "study time := 00:00:00\n";
    output_h33 << "imagedata byte order := LITTLEENDIAN\n;\n";
    output_h33 << "number of energy windows := 1\n;\n";

    output_h33 << "energy window [1] :=\n";
    output_h33 << "energy window lower level [1] :=\n;";
    output_h33 << "energy window upper level [1] :=\n";
    output_h33 << "flood corrected := N\n";
    output_h33 << "decay corrected := N\n;\n";
    output_h33 << "!SPECT STUDY (general) :=\n";
    output_h33 << "number of detector heads := 1\n;\n";

    output_h33 << "!number of images/energy window := " << image.Dimensions()[2] << std::endl;
    output_h33 << "!process status := Reconstructed\n";
    output_h33 << "!matrix size [1] := " << image.Dimensions()[0] << std::endl;
    output_h33 << "!matrix size [2] := " << image.Dimensions()[1] << std::endl;
    output_h33 << "!number format := " << number_format << std::endl;
    output_h33 << "!number of bytes per pixel := " << number_of_bytes_per_pixel << std::endl;
    output_h33 << "scaling factor (mm/pixel) [1] := +" << image.Spacing()[0] << std::endl;
    output_h33 << "scaling factor (mm/pixel) [2] := +" << image.Spacing()[1] << std::endl;
    output_h33 << "!number of projections := " << image.Dimensions()[2] << std::endl;
    output_h33 << "!extent of rotation :=\n";
    output_h33 << "!time per projection (sec) := 0\n";
    output_h33 << "study duration (sec) := 0\n";
    output_h33 << "!maximum pixel count := " << image.NumVoxels() << std::endl;
    output_h33 << "patient orientation := Unknown\n";
    output_h33 << "patient rotation := Unknown\n;\n";

    output_h33 << "!SPECT STUDY (reconstructed data) :=\n";
    output_h33 << "method of reconstruction := Unknown\n";
    output_h33 << "!number of slices := " << image.Dimensions()[2] << std::endl;
    output_h33 << "number of reference frame := 0\n";
    output_h33 << "slice orientation := Transverse\n";
    output_h33 << "slice thickness (pixels) := +" << image.Spacing()[2] << std::endl;
    output_h33 << "centre-centre slice separation (pixels) := +" << image.Spacing()[2] << std::endl;
    output_h33 << "filter name := Unknown\n";
    output_h33 << "filter parameters := Cutoff\n";
    output_h33 << "method of attenuation correction := measured\n";
    output_h33 << "scatter corrected := N\n";
    output_h33 << "oblique reconstruction := N\n";
    output_h33 << "!END OF INTERFILE :=\n";

    // Close header file.
    output_h33.close();

}


template<typename VOXELTYPE>
void H33IO::ReadI33Data(const std::string &i33_filename, VoxelImage<VOXELTYPE> &image)
{
    // Load binary file.
    std::ifstream binary_file(i33_filename.c_str(), std::ios::in | std::ios::binary);

    // Check if binary file is opened.
    if (!binary_file.is_open()) {
        std::string error = "[IMP ERROR] Cannot open binary interfile image (.i33).";
        throw std::invalid_argument(error.c_str());
    }

    //Check file extension.
    std::string extension = i33_filename.substr(i33_filename.find_last_of("."));
    if (extension != ".i33") {
        std::string error = "[IMP ERROR] Cannot read binary interfile image of unknown format. Expected format: (.i33)";
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


// template<typename VOXELTYPE>
// void H33IO::ReadI33GzData(const std::string &i33_gz_filename, VoxelImage<VOXELTYPE> &image)
// {
//     //Check file extension.
//     std::string extension = i33_gz_filename.substr(i33_gz_filename.find_last_of("."));
//     if (extension != ".gz") {
//         std::string error = "[IMP ERROR] Cannot read compressed binary interfile (.i33) of unknown compression. Expected compression: (.gz)";
//         throw std::invalid_argument(error.c_str());
//     }

//     // Check if image is initialized.
//     if (!image.IsInitialized()) {
//         std::string error = "[IMP ERROR] The voxelized image is not initialized. Set dimensions, number of voxels and memory size.";
//         throw std::runtime_error(error.c_str());
//     }

//     // Clear the image data.
//     image.EditData().clear();
//     image.EditData().reserve(image.NumVoxels());

//     // Compressed file reading scope.
//     {
//         // Read compressed file.
//         std::ifstream file_in(i33_gz_filename, std::ios::in | std::ios_base::binary);
//         boost::iostreams::filtering_istream input;
//         input.push(boost::iostreams::gzip_decompressor());
//         input.push(file_in);

//         // Get character pointer to the image data vector.
//         char *data_ptr = reinterpret_cast<char*>(image.EditData().data());

//         // Assign data to the image data vector byte per byte.
//         for (int i = 0; i != image.MemorySize(); ++i) {
//             input.read(&data_ptr[i], 1);
//         }

//         // Close compressed binary file.
//         file_in.close();

//         // Release character pointer.
//         data_ptr = nullptr;

//     } // End of compressed file reading scope.

// }


template<typename VOXELTYPE>
void H33IO::WriteI33Data(VoxelImage<VOXELTYPE> &image, const std::string &i33_filename)
{
    // Create the binary file.
    std::ofstream output_i33(i33_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

    if (!output_i33.is_open()) {
        std::string error = "[IMP ERROR] Could not open binary file to save voxelized image in interfile format (.h33).";
        throw std::invalid_argument(error.c_str());
    }

    // Write binary data.
    output_i33.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

    // Close binary file.
    output_i33.close();

}


// template<typename VOXELTYPE>
// void H33IO::WriteI33GzData(VoxelImage<VOXELTYPE> &image, const std::string &i33_gz_filename)
// {
//     // Create the compressed binary file.
//     std::ofstream output_i33_gz(i33_gz_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

//     // Check if compressed binary file can be opened.
//     if (!output_i33_gz.is_open()) {
//         std::string error = "[IMP ERROR] Could not open binary file to save voxelized image in interfile format (.h33).";
//         throw std::invalid_argument(error.c_str());
//     }

//     // Convert image binary data to sting stream.
//     std::stringstream uncompressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);
//     uncompressed_data.write(reinterpret_cast<const char*>(&image.Data()[0]), image.MemorySize());

//     // Sting stream to store the compressed binary data.
//     std::stringstream compressed_data(std::stringstream::in | std::stringstream::out | std::stringstream::binary);

//     // Set the image's compressed state ON.
//     image.SetCompression(ImageCompression::ON);

//     // Compressed file writing scope.
//     {
//         // Create compressed output buffer.
//         boost::iostreams::filtering_streambuf<boost::iostreams::input> out_buffer;
//         out_buffer.push(boost::iostreams::counter());
//         out_buffer.push(boost::iostreams::gzip_compressor());
//         out_buffer.push(uncompressed_data);

//         // Copy compressed output buffer in compressed binary file.
//         boost::iostreams::copy(out_buffer, compressed_data);

//         // Flush the data from the compressor.
//         boost::iostreams::close(out_buffer);

//         // Set the number of the compressed bytes of the image.
//         image.SetCompressedBytes(compressed_data.str().length());

//     }

//     // Write the compressed image binary data.
//     output_i33_gz.write(compressed_data.str().c_str(), compressed_data.str().length());

//     // Close compressed binary file.
//     output_i33_gz.close();

// }


template <typename DATATYPE>
void H33IO::ExtractAttributeFromHeader(const std::string &line, DATATYPE &attribute)
{
    // Convert line string to stringstream.
    std::stringstream ss(line);

    // First part of header's line.
    std::string token = "";

    // Get the attribute by processing the stringstream.
    if (!(ss >> token >> attribute)) {
        std::string error = "[IMP ERROR] Can not extract attribute from the Interfile image header (.h33) from line: " + line;
        throw std::runtime_error(error.c_str());
    }

}


} //end of namespace IMP


#endif //IMP_ENGINE_IMAGE_IO_H33_IO_TPP_
