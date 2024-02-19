#include "voxel_image_test.hpp"

namespace IMP_test {


// Check default constructor.
TEST_F(VoxelImageTest, DefaultCtor) {

    // Check size of image data.
    EXPECT_EQ(0, this->uint_image_.Data().size());
    EXPECT_EQ(0, this->int_image_.Data().size());
    EXPECT_EQ(0, this->float_image_.Data().size());
    EXPECT_EQ(0, this->double_image_.Data().size());

    // Check initialization of image dimensions.
    EXPECT_EQ(-1, this->uint_image_.Dimensions().At(0));
    EXPECT_EQ(-1, this->uint_image_.Dimensions().At(1));
    EXPECT_EQ(-1, this->uint_image_.Dimensions().At(2));
    //
    EXPECT_EQ(-1, this->int_image_.Dimensions().At(0));
    EXPECT_EQ(-1, this->int_image_.Dimensions().At(1));
    EXPECT_EQ(-1, this->int_image_.Dimensions().At(2));
    //
    EXPECT_EQ(-1, this->float_image_.Dimensions().At(0));
    EXPECT_EQ(-1, this->float_image_.Dimensions().At(1));
    EXPECT_EQ(-1, this->float_image_.Dimensions().At(2));
    //
    EXPECT_EQ(-1, this->double_image_.Dimensions().At(0));
    EXPECT_EQ(-1, this->double_image_.Dimensions().At(1));
    EXPECT_EQ(-1, this->double_image_.Dimensions().At(2));

    // Check initialization of image spacing.
    EXPECT_DOUBLE_EQ(-1., this->uint_image_.Spacing().At(0));
    EXPECT_DOUBLE_EQ(-1., this->uint_image_.Spacing().At(1));
    EXPECT_DOUBLE_EQ(-1., this->uint_image_.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(-1., this->int_image_.Spacing().At(0));
    EXPECT_DOUBLE_EQ(-1., this->int_image_.Spacing().At(1));
    EXPECT_DOUBLE_EQ(-1., this->int_image_.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(-1., this->float_image_.Spacing().At(0));
    EXPECT_DOUBLE_EQ(-1., this->float_image_.Spacing().At(1));
    EXPECT_DOUBLE_EQ(-1., this->float_image_.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(-1., this->double_image_.Spacing().At(0));
    EXPECT_DOUBLE_EQ(-1., this->double_image_.Spacing().At(1));
    EXPECT_DOUBLE_EQ(-1., this->double_image_.Spacing().At(2));

    // Check initialization of image Origin.
    EXPECT_DOUBLE_EQ(0., this->uint_image_.Origin().At(0));
    EXPECT_DOUBLE_EQ(0., this->uint_image_.Origin().At(1));
    EXPECT_DOUBLE_EQ(0., this->uint_image_.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(0., this->int_image_.Origin().At(0));
    EXPECT_DOUBLE_EQ(0., this->int_image_.Origin().At(1));
    EXPECT_DOUBLE_EQ(0., this->int_image_.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(0., this->float_image_.Origin().At(0));
    EXPECT_DOUBLE_EQ(0., this->float_image_.Origin().At(1));
    EXPECT_DOUBLE_EQ(0., this->float_image_.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(0., this->double_image_.Origin().At(0));
    EXPECT_DOUBLE_EQ(0., this->double_image_.Origin().At(1));
    EXPECT_DOUBLE_EQ(0., this->double_image_.Origin().At(2));

    // Check initialization of image number of voxels.
    EXPECT_EQ(-1, this->uint_image_.NumVoxels());
    EXPECT_EQ(-1, this->int_image_.NumVoxels());
    EXPECT_EQ(-1, this->float_image_.NumVoxels());
    EXPECT_EQ(-1, this->double_image_.NumVoxels());

    // Check initialization of image memory size.
    EXPECT_EQ(-1, this->uint_image_.MemorySize());
    EXPECT_EQ(-1, this->int_image_.MemorySize());
    EXPECT_EQ(-1, this->float_image_.MemorySize());
    EXPECT_EQ(-1, this->double_image_.MemorySize());

    // Check initialization of image data type.
    EXPECT_EQ(IMP::ImageDataType::no_type, this->uint_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::no_type, this->int_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::no_type, this->float_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::no_type, this->double_image_.DataType());

    // Check initialization of image compression state.
    EXPECT_EQ(IMP::ImageCompression::OFF, this->uint_image_.Compression());
    EXPECT_EQ(IMP::ImageCompression::OFF, this->int_image_.Compression());
    EXPECT_EQ(IMP::ImageCompression::OFF, this->float_image_.Compression());
    EXPECT_EQ(IMP::ImageCompression::OFF, this->double_image_.Compression());

    // Check initialization of image compression bytes.
    EXPECT_EQ(-1, this->uint_image_.CompressedBytes());
    EXPECT_EQ(-1, this->int_image_.CompressedBytes());
    EXPECT_EQ(-1, this->float_image_.CompressedBytes());
    EXPECT_EQ(-1, this->double_image_.CompressedBytes());

}


// Check copy constructor.
TEST_F(VoxelImageTest, CopyCtor) {

    // Load test image datasets.
    IMP::VoxelImage<unsigned> test_img_uint;
    test_img_uint.LoadFrom("../test_data/test_img_uint.mhd");
    IMP::VoxelImage<int> test_img_int;
    test_img_int.LoadFrom("../test_data/test_img_int.mhd");
    IMP::VoxelImage<float> test_img_float;
    test_img_float.LoadFrom("../test_data/test_img_float.mhd");
    IMP::VoxelImage<double> test_img_double;
    test_img_double.LoadFrom("../test_data/test_img_double.mhd");

    // Create new images as copies of the loaded ones.
    IMP::VoxelImage<unsigned> img_uint(test_img_uint);
    IMP::VoxelImage<int> img_int(test_img_int);
    IMP::VoxelImage<float> img_float(test_img_float);
    IMP::VoxelImage<double> img_double(test_img_double);

    // Check image dimensions consistency.
    EXPECT_EQ(test_img_uint.Dimensions().At(0), img_uint.Dimensions().At(0));
    EXPECT_EQ(test_img_uint.Dimensions().At(1), img_uint.Dimensions().At(1));
    EXPECT_EQ(test_img_uint.Dimensions().At(2), img_uint.Dimensions().At(2));
    //
    EXPECT_EQ(test_img_int.Dimensions().At(0), img_int.Dimensions().At(0));
    EXPECT_EQ(test_img_int.Dimensions().At(1), img_int.Dimensions().At(1));
    EXPECT_EQ(test_img_uint.Dimensions().At(2), img_int.Dimensions().At(2));
    //
    EXPECT_EQ(test_img_float.Dimensions().At(0), img_float.Dimensions().At(0));
    EXPECT_EQ(test_img_float.Dimensions().At(1), img_float.Dimensions().At(1));
    EXPECT_EQ(test_img_float.Dimensions().At(2), img_float.Dimensions().At(2));
    //
    EXPECT_EQ(test_img_double.Dimensions().At(0), img_double.Dimensions().At(0));
    EXPECT_EQ(test_img_double.Dimensions().At(1), img_double.Dimensions().At(1));
    EXPECT_EQ(test_img_double.Dimensions().At(2), img_double.Dimensions().At(2));

    // Check image spacing consistency.
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(0), img_uint.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(1), img_uint.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(2), img_uint.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_int.Spacing().At(0), img_int.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_int.Spacing().At(1), img_int.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(2), img_int.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_float.Spacing().At(0), img_float.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_float.Spacing().At(1), img_float.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_float.Spacing().At(2), img_float.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_double.Spacing().At(0), img_double.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_double.Spacing().At(1), img_double.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_double.Spacing().At(2), img_double.Spacing().At(2));

    // Check image Origin consistency.
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(0), img_uint.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(1), img_uint.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(2), img_uint.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_int.Origin().At(0), img_int.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_int.Origin().At(1), img_int.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(2), img_int.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_float.Origin().At(0), img_float.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_float.Origin().At(1), img_float.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_float.Origin().At(2), img_float.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_double.Origin().At(0), img_double.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_double.Origin().At(1), img_double.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_double.Origin().At(2), img_double.Origin().At(2));

    // Check image number of voxels consistency.
    EXPECT_EQ(test_img_uint.NumVoxels(), img_uint.NumVoxels());
    EXPECT_EQ(test_img_int.NumVoxels(), img_int.NumVoxels());
    EXPECT_EQ(test_img_float.NumVoxels(), img_float.NumVoxels());
    EXPECT_EQ(test_img_double.NumVoxels(), img_double.NumVoxels());

    // Check image memory size consistency.
    EXPECT_EQ(test_img_uint.MemorySize(), img_uint.MemorySize());
    EXPECT_EQ(test_img_int.MemorySize(), img_int.MemorySize());
    EXPECT_EQ(test_img_float.MemorySize(), img_float.MemorySize());
    EXPECT_EQ(test_img_double.MemorySize(), img_double.MemorySize());

    // Check image data type consistency.
    EXPECT_EQ(test_img_uint.DataType(), img_uint.DataType());
    EXPECT_EQ(test_img_int.DataType(), img_int.DataType());
    EXPECT_EQ(test_img_float.DataType(), img_float.DataType());
    EXPECT_EQ(test_img_double.DataType(), img_double.DataType());

    // Check image compression state consistency.
    EXPECT_EQ(test_img_uint.Compression(), img_uint.Compression());
    EXPECT_EQ(test_img_int.Compression(), img_int.Compression());
    EXPECT_EQ(test_img_float.Compression(), img_float.Compression());
    EXPECT_EQ(test_img_double.Compression(), img_double.Compression());

    // Check image compression bytes consistency.
    EXPECT_EQ(test_img_uint.CompressedBytes(), img_uint.CompressedBytes());
    EXPECT_EQ(test_img_int.CompressedBytes(), img_int.CompressedBytes());
    EXPECT_EQ(test_img_float.CompressedBytes(), img_float.CompressedBytes());
    EXPECT_EQ(test_img_double.CompressedBytes(), img_double.CompressedBytes());

    // Check data consistency.
    for (int i = 0; i < img_uint.NumVoxels(); ++i) {
        EXPECT_EQ(test_img_uint.Data()[i], img_uint.Data()[i]);
        EXPECT_EQ(test_img_int.Data()[i], img_int.Data()[i]);
        EXPECT_FLOAT_EQ(test_img_float.Data()[i], img_float.Data()[i]);
        EXPECT_DOUBLE_EQ(test_img_double.Data()[i], img_double.Data()[i]);
    }
}


// Check move constructor.
TEST_F(VoxelImageTest, MoveCtor) {

    // Load test image datasets.
    IMP::VoxelImage<unsigned> test_img_uint;
    test_img_uint.LoadFrom("../test_data/test_img_uint.mhd");
    IMP::VoxelImage<int> test_img_int;
    test_img_int.LoadFrom("../test_data/test_img_int.mhd");
    IMP::VoxelImage<float> test_img_float;
    test_img_float.LoadFrom("../test_data/test_img_float.mhd");
    IMP::VoxelImage<double> test_img_double;
    test_img_double.LoadFrom("../test_data/test_img_double.mhd");

    // Create copy images of the loaded ones.
    IMP::VoxelImage<unsigned> cp_test_img_uint(test_img_uint);
    IMP::VoxelImage<int> cp_test_img_int(test_img_int);
    IMP::VoxelImage<float> cp_test_img_float(test_img_float);
    IMP::VoxelImage<double> cp_test_img_double(test_img_double);

    // Create new images by moving the copy images.
    IMP::VoxelImage<unsigned> img_uint(std::move(cp_test_img_uint));
    IMP::VoxelImage<int> img_int(std::move(cp_test_img_int));
    IMP::VoxelImage<float> img_float(std::move(cp_test_img_float));
    IMP::VoxelImage<double> img_double(std::move(cp_test_img_double));

    // Check image dimensions consistency.
    EXPECT_EQ(test_img_uint.Dimensions().At(0), img_uint.Dimensions().At(0));
    EXPECT_EQ(test_img_uint.Dimensions().At(1), img_uint.Dimensions().At(1));
    EXPECT_EQ(test_img_uint.Dimensions().At(2), img_uint.Dimensions().At(2));
    //
    EXPECT_EQ(test_img_int.Dimensions().At(0), img_int.Dimensions().At(0));
    EXPECT_EQ(test_img_int.Dimensions().At(1), img_int.Dimensions().At(1));
    EXPECT_EQ(test_img_uint.Dimensions().At(2), img_int.Dimensions().At(2));
    //
    EXPECT_EQ(test_img_float.Dimensions().At(0), img_float.Dimensions().At(0));
    EXPECT_EQ(test_img_float.Dimensions().At(1), img_float.Dimensions().At(1));
    EXPECT_EQ(test_img_float.Dimensions().At(2), img_float.Dimensions().At(2));
    //
    EXPECT_EQ(test_img_double.Dimensions().At(0), img_double.Dimensions().At(0));
    EXPECT_EQ(test_img_double.Dimensions().At(1), img_double.Dimensions().At(1));
    EXPECT_EQ(test_img_double.Dimensions().At(2), img_double.Dimensions().At(2));

    // Check image spacing consistency.
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(0), img_uint.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(1), img_uint.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(2), img_uint.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_int.Spacing().At(0), img_int.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_int.Spacing().At(1), img_int.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Spacing().At(2), img_int.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_float.Spacing().At(0), img_float.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_float.Spacing().At(1), img_float.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_float.Spacing().At(2), img_float.Spacing().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_double.Spacing().At(0), img_double.Spacing().At(0));
    EXPECT_DOUBLE_EQ(test_img_double.Spacing().At(1), img_double.Spacing().At(1));
    EXPECT_DOUBLE_EQ(test_img_double.Spacing().At(2), img_double.Spacing().At(2));

    // Check image Origin consistency.
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(0), img_uint.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(1), img_uint.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(2), img_uint.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_int.Origin().At(0), img_int.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_int.Origin().At(1), img_int.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_uint.Origin().At(2), img_int.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_float.Origin().At(0), img_float.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_float.Origin().At(1), img_float.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_float.Origin().At(2), img_float.Origin().At(2));
    //
    EXPECT_DOUBLE_EQ(test_img_double.Origin().At(0), img_double.Origin().At(0));
    EXPECT_DOUBLE_EQ(test_img_double.Origin().At(1), img_double.Origin().At(1));
    EXPECT_DOUBLE_EQ(test_img_double.Origin().At(2), img_double.Origin().At(2));

    // Check image number of voxels consistency.
    EXPECT_EQ(test_img_uint.NumVoxels(), img_uint.NumVoxels());
    EXPECT_EQ(test_img_int.NumVoxels(), img_int.NumVoxels());
    EXPECT_EQ(test_img_float.NumVoxels(), img_float.NumVoxels());
    EXPECT_EQ(test_img_double.NumVoxels(), img_double.NumVoxels());

    // Check image memory size consistency.
    EXPECT_EQ(test_img_uint.MemorySize(), img_uint.MemorySize());
    EXPECT_EQ(test_img_int.MemorySize(), img_int.MemorySize());
    EXPECT_EQ(test_img_float.MemorySize(), img_float.MemorySize());
    EXPECT_EQ(test_img_double.MemorySize(), img_double.MemorySize());

    // Check image data type consistency.
    EXPECT_EQ(test_img_uint.DataType(), img_uint.DataType());
    EXPECT_EQ(test_img_int.DataType(), img_int.DataType());
    EXPECT_EQ(test_img_float.DataType(), img_float.DataType());
    EXPECT_EQ(test_img_double.DataType(), img_double.DataType());

    // Check image compression state consistency.
    EXPECT_EQ(test_img_uint.Compression(), img_uint.Compression());
    EXPECT_EQ(test_img_int.Compression(), img_int.Compression());
    EXPECT_EQ(test_img_float.Compression(), img_float.Compression());
    EXPECT_EQ(test_img_double.Compression(), img_double.Compression());

    // Check image compression bytes consistency.
    EXPECT_EQ(test_img_uint.CompressedBytes(), img_uint.CompressedBytes());
    EXPECT_EQ(test_img_int.CompressedBytes(), img_int.CompressedBytes());
    EXPECT_EQ(test_img_float.CompressedBytes(), img_float.CompressedBytes());
    EXPECT_EQ(test_img_double.CompressedBytes(), img_double.CompressedBytes());

    // Check data consistency.
    for (int i = 0; i < img_uint.NumVoxels(); ++i) {
        EXPECT_EQ(test_img_uint.Data()[i], img_uint.Data()[i]);
        EXPECT_EQ(test_img_int.Data()[i], img_int.Data()[i]);
        EXPECT_FLOAT_EQ(test_img_float.Data()[i], img_float.Data()[i]);
        EXPECT_DOUBLE_EQ(test_img_double.Data()[i], img_double.Data()[i]);
    }

}


// Check InitEmptyImage function.
TEST_F(VoxelImageTest, InitEmptyImage) {

    // Set image dimensions, spacing, and Origin.
    IMP::Vec<3,int> dims({10,10,10});
    IMP::Vec<3,double> spacing({0.4,0.4,2});
    IMP::Vec<3,double> Origin({-0.5,3,15.2});

    // Create empty images.
    IMP::VoxelImage<unsigned> img_uint; img_uint.InitEmptyImage(dims, spacing, Origin);
    IMP::VoxelImage<int> img_int; img_int.InitEmptyImage(dims, spacing, Origin);
    IMP::VoxelImage<float> img_float; img_float.InitEmptyImage(dims, spacing, Origin);
    IMP::VoxelImage<double> img_double; img_double.InitEmptyImage(dims, spacing, Origin);

    // Check the initialization of the images dimensions.
    EXPECT_EQ(10, img_uint.Dimensions()[0]); EXPECT_EQ(10, img_uint.Dimensions()[1]); EXPECT_EQ(10, img_uint.Dimensions()[2]);
    EXPECT_EQ(10, img_int.Dimensions()[0]); EXPECT_EQ(10, img_int.Dimensions()[1]); EXPECT_EQ(10, img_int.Dimensions()[2]);
    EXPECT_EQ(10, img_float.Dimensions()[0]); EXPECT_EQ(10, img_float.Dimensions()[1]); EXPECT_EQ(10, img_float.Dimensions()[2]);
    EXPECT_EQ(10, img_double.Dimensions()[0]); EXPECT_EQ(10, img_double.Dimensions()[1]); EXPECT_EQ(10, img_double.Dimensions()[2]);

    // Check the initialization of the images spacing.
    EXPECT_DOUBLE_EQ(0.4, img_uint.Spacing()[0]); EXPECT_EQ(0.4, img_uint.Spacing()[1]); EXPECT_EQ(2., img_uint.Spacing()[2]);
    EXPECT_DOUBLE_EQ(0.4, img_int.Spacing()[0]); EXPECT_EQ(0.4, img_int.Spacing()[1]); EXPECT_EQ(2., img_int.Spacing()[2]);
    EXPECT_DOUBLE_EQ(0.4, img_float.Spacing()[0]); EXPECT_EQ(0.4, img_float.Spacing()[1]); EXPECT_EQ(2., img_float.Spacing()[2]);
    EXPECT_DOUBLE_EQ(0.4, img_double.Spacing()[0]); EXPECT_EQ(0.4, img_double.Spacing()[1]); EXPECT_EQ(2., img_double.Spacing()[2]);

    // Check the initialization of the images Origin.
    EXPECT_DOUBLE_EQ(-0.5, img_uint.Origin()[0]); EXPECT_EQ(3, img_uint.Origin()[1]); EXPECT_EQ(15.2, img_uint.Origin()[2]);
    EXPECT_DOUBLE_EQ(-0.5, img_int.Origin()[0]); EXPECT_EQ(3, img_int.Origin()[1]); EXPECT_EQ(15.2, img_int.Origin()[2]);
    EXPECT_DOUBLE_EQ(-0.5, img_float.Origin()[0]); EXPECT_EQ(3, img_float.Origin()[1]); EXPECT_EQ(15.2, img_float.Origin()[2]);
    EXPECT_DOUBLE_EQ(-0.5, img_double.Origin()[0]); EXPECT_EQ(3, img_double.Origin()[1]); EXPECT_EQ(15.2, img_double.Origin()[2]);

    // Check the initialization of the images number of voxels.
    int num_voxels = dims[0]*dims[1]*dims[2];
    EXPECT_EQ(num_voxels, img_uint.NumVoxels()); EXPECT_EQ(num_voxels, img_int.NumVoxels());
    EXPECT_EQ(num_voxels, img_float.NumVoxels()); EXPECT_EQ(num_voxels, img_double.NumVoxels());

    // Check the initialization of the images memory size.
    EXPECT_EQ(num_voxels*sizeof(unsigned), img_uint.MemorySize()); EXPECT_EQ(num_voxels*sizeof(int), img_int.MemorySize());
    EXPECT_EQ(num_voxels*sizeof(float), img_float.MemorySize()); EXPECT_EQ(num_voxels*sizeof(double), img_double.MemorySize());

    // Check the initialization of the images type.
    EXPECT_EQ(IMP::ImageDataType::no_type, img_uint.DataType()); EXPECT_EQ(IMP::ImageDataType::no_type, img_int.DataType());
    EXPECT_EQ(IMP::ImageDataType::no_type, img_float.DataType()); EXPECT_EQ(IMP::ImageDataType::no_type, img_double.DataType());

    // Check the initialization of the compression status.
    EXPECT_EQ(IMP::ImageCompression::OFF, img_uint.Compression()); EXPECT_EQ(IMP::ImageCompression::OFF, img_int.Compression());
    EXPECT_EQ(IMP::ImageCompression::OFF, img_float.Compression()); EXPECT_EQ(IMP::ImageCompression::OFF, img_double.Compression());

    // Check the initialization of the compression bytes number.
    EXPECT_EQ(-1, img_uint.CompressedBytes()); EXPECT_EQ(-1, img_int.CompressedBytes());
    EXPECT_EQ(-1, img_float.CompressedBytes()); EXPECT_EQ(-1, img_double.CompressedBytes());

    // Check the initialization of the images data.
    for (int i = 0; i < img_uint.NumVoxels(); ++i) {
        EXPECT_EQ(0, img_uint.Data()[i]);
        EXPECT_EQ(0, img_int.Data()[i]);
        EXPECT_FLOAT_EQ(0.f, img_float.Data()[i]);
        EXPECT_DOUBLE_EQ(0., img_double.Data()[i]);
    }

}


// Check SetDimensions function.
TEST_F(VoxelImageTest, SetDimensions) {

    // Set meaningfull dimensions for uint_image.
    this->uint_image_.SetDimensions(10,10,5);

    // Check the assigned dimensions.
    EXPECT_EQ(10, this->uint_image_.Dimensions()[0]);
    EXPECT_EQ(10, this->uint_image_.Dimensions()[1]);
    EXPECT_EQ(5, this->uint_image_.Dimensions()[2]);

    // Set meaningless dimensions and check for invalid argument exception.
    EXPECT_THROW(this->uint_image_.SetDimensions(-10,10,5), std::invalid_argument);
    EXPECT_THROW(this->int_image_.SetDimensions(10,-10,5), std::invalid_argument);
    EXPECT_THROW(this->float_image_.SetDimensions(10,10,-5), std::invalid_argument);
    EXPECT_THROW(this->double_image_.SetDimensions(-10,-10,-5), std::invalid_argument);

}


// Check SetSpacing function.
TEST_F(VoxelImageTest, SetSpacing) {

    // Set meaningfull spacing for uint_image.
    this->uint_image_.SetSpacing(0.15,1.45,150);

    // Check the assigned spacing.
    EXPECT_DOUBLE_EQ(0.15, this->uint_image_.Spacing()[0]);
    EXPECT_DOUBLE_EQ(1.45, this->uint_image_.Spacing()[1]);
    EXPECT_DOUBLE_EQ(150, this->uint_image_.Spacing()[2]);

    // Set meaningless spacing and check for invalid argument exception.
    EXPECT_THROW(this->uint_image_.SetSpacing(-10,10,5), std::invalid_argument);
    EXPECT_THROW(this->int_image_.SetSpacing(-0.000001,-0,5), std::invalid_argument);
    EXPECT_THROW(this->float_image_.SetSpacing(-0.0001,10,-5), std::invalid_argument);

    // Check that doesn't throw at the voxel limit.
    EXPECT_NO_THROW(this->double_image_.SetSpacing(0.000001,0.000001,0.000001));

}


// Check SetOrigin function.
TEST_F(VoxelImageTest, SetOrigin) {

    // Set meaningfull Origin for uint_image.
    this->uint_image_.SetOrigin(0.15,1.45,150);

    // Check the assigned Origin.
    EXPECT_DOUBLE_EQ(0.15, this->uint_image_.Origin()[0]);
    EXPECT_DOUBLE_EQ(1.45, this->uint_image_.Origin()[1]);
    EXPECT_DOUBLE_EQ(150, this->uint_image_.Origin()[2]);

    // Set meaningfull Origin for int_image.
    this->int_image_.SetOrigin(-0.15,-1.45,-150);

    // Check the assigned Origin.
    EXPECT_DOUBLE_EQ(-0.15, this->int_image_.Origin()[0]);
    EXPECT_DOUBLE_EQ(-1.45, this->int_image_.Origin()[1]);
    EXPECT_DOUBLE_EQ(-150, this->int_image_.Origin()[2]);

}


// Check SetNumVoxels function.
TEST_F(VoxelImageTest, SetNumVoxels) {

    // Set meaningfull number of voxels.
    this->uint_image_.SetNumVoxels(0);
    this->int_image_.SetNumVoxels(150);
    this->float_image_.SetNumVoxels(-0);

    // Check the assigned number of voxels.
    EXPECT_EQ(0, this->uint_image_.NumVoxels());
    EXPECT_EQ(150, this->int_image_.NumVoxels());
    EXPECT_EQ(0, this->float_image_.NumVoxels());

    // Set meaningless number of voxels and check for invalid argument exception.
    EXPECT_THROW(this->int_image_.SetNumVoxels(-12), std::invalid_argument);

}


// Check SetMemorySize function.
TEST_F(VoxelImageTest, SetMemorySize) {

    // Set meaningfull memory size.
    this->uint_image_.SetMemorySize(0);
    this->int_image_.SetMemorySize(150);
    this->float_image_.SetMemorySize(-0);

    // Check the assigned memory size.
    EXPECT_EQ(0, this->uint_image_.MemorySize());
    EXPECT_EQ(150, this->int_image_.MemorySize());
    EXPECT_EQ(0, this->float_image_.MemorySize());

    // Set meaningless number of voxels and check for invalid argument exception.
    EXPECT_THROW(this->int_image_.SetMemorySize(-12), std::invalid_argument);

}


// Check SetDataType function.
TEST_F(VoxelImageTest, SetDataType) {

    // Set image data types.
    this->uint_image_.SetDataType(IMP::ImageDataType::no_type);
    this->int_image_.SetDataType(IMP::ImageDataType::char_type);
    this->float_image_.SetDataType(IMP::ImageDataType::uchar_type);
    this->double_image_.SetDataType(IMP::ImageDataType::short_int_type);

    // Check the assigned data types.
    EXPECT_EQ(IMP::ImageDataType::no_type, this->uint_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::char_type, this->int_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::uchar_type, this->float_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::short_int_type, this->double_image_.DataType());

    // Set more image data types.
    this->uint_image_.SetDataType(IMP::ImageDataType::ushort_int_type);
    this->int_image_.SetDataType(IMP::ImageDataType::int_type);
    this->float_image_.SetDataType(IMP::ImageDataType::uint_type);
    this->double_image_.SetDataType(IMP::ImageDataType::long_int_type);

    // Check the assigned data types.
    EXPECT_EQ(IMP::ImageDataType::ushort_int_type, this->uint_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::int_type, this->int_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::uint_type, this->float_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::long_int_type, this->double_image_.DataType());

    // Set more image data types.
    this->uint_image_.SetDataType(IMP::ImageDataType::ulong_int_type);
    this->int_image_.SetDataType(IMP::ImageDataType::float_type);
    this->float_image_.SetDataType(IMP::ImageDataType::double_type);
    this->double_image_.SetDataType(IMP::ImageDataType::bool_type);

    // Check the assigned data types.
    EXPECT_EQ(IMP::ImageDataType::ulong_int_type, this->uint_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::float_type, this->int_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::double_type, this->float_image_.DataType());
    EXPECT_EQ(IMP::ImageDataType::bool_type, this->double_image_.DataType());

}


// Check SetCompression function.
TEST_F(VoxelImageTest, SetCompression) {

    // Set image compression types.
    this->uint_image_.SetCompression(IMP::ImageCompression::ON);
    this->int_image_.SetCompression(IMP::ImageCompression::OFF);

    // Check the assigned compression types.
    EXPECT_EQ(IMP::ImageCompression::ON, this->uint_image_.Compression());
    EXPECT_EQ(IMP::ImageCompression::OFF, this->int_image_.Compression());

}


// Check SetCompression function.
TEST_F(VoxelImageTest, SetCompressedBytes) {

    // Set image compression types.
    this->uint_image_.SetCompressedBytes(5);

    // Check the assigned compression types.
    EXPECT_EQ(5, this->uint_image_.CompressedBytes());
    EXPECT_THROW(this->int_image_.SetCompressedBytes(-5), std::invalid_argument);

}


} // End of namespace IMP_test.
