#include "vec_test.hpp"

namespace IMP_test {


// Check if data of vector is properly initialized.
TEST_F(VecTest, DefaultCtor) {

    // Check size of vectors' data.
    EXPECT_EQ(1, this->vec_u1_.Dim());
    EXPECT_EQ(2, this->vec_i2_.Dim());
    EXPECT_EQ(3, this->vec_f3_.Dim());
    EXPECT_EQ(4, this->vec_d4_.Dim());

    // Check vectors' data zero.
    for (const auto &val_u1 : this->vec_u1_) { EXPECT_EQ(0, val_u1); }
    for (const auto &val_i2 : this->vec_i2_) { EXPECT_EQ(0, val_i2); }
    for (const auto &val_f3 : this->vec_f3_) { EXPECT_EQ(0, val_f3); }
    for (const auto &val_d4 : this->vec_d4_) { EXPECT_EQ(0, val_d4); }

}


// Check if dimensions of vector corresponds to data size after initialization.
TEST_F(VecTest, DimensionInit) {

    // Check size of vectors' data with vectors' dim.
    EXPECT_EQ(1, this->vec_u1_.Dim());
    EXPECT_EQ(2, this->vec_i2_.Dim());
    EXPECT_EQ(3, this->vec_f3_.Dim());
    EXPECT_EQ(4, this->vec_d4_.Dim());

}


// Check copy constructor.
TEST_F(VecTest, CopyCtor) {

    // Create a vector.
    IMP::Vec<3, double> vec1; vec1.Set({12., 15., -30.});

    // Create new vector by copying previous.
    IMP::Vec<3, double> vec2(vec1);

    // Check.
    EXPECT_TRUE(vec1 == vec2);

}


// Check move constructor.
TEST_F(VecTest, MoveCtor) {

    // Create a vector.
    IMP::Vec<3, double> vec1; vec1.Set({12., 15., -30.});

    // Create new vector by moving previous.
    IMP::Vec<3, double> vec2(std::move(vec1));

    // Check.
    EXPECT_DOUBLE_EQ(12., vec2.At(0));
    EXPECT_DOUBLE_EQ(15., vec2.At(1));
    EXPECT_DOUBLE_EQ(-30., vec2.At(2));
    EXPECT_EQ(0, std::distance(vec1.begin(), vec1.end()) );

}


// Check constructor from std::iniatilizer_list. More extensive test in the corresponding Set() function.
TEST_F(VecTest, InitListCtor) {

    // Set to random values.
    IMP::Vec<3,double> vec({5., -3., 10.});

    // Check.
    EXPECT_DOUBLE_EQ(5., vec.At(0));
    EXPECT_DOUBLE_EQ(-3., vec.At(1));
    EXPECT_DOUBLE_EQ(10., vec.At(2));

}

// Check constructor from std::vector. More extensive test in the corresponding Set() function.
TEST_F(VecTest, InitFromStdVector) {

    // Set to random values.
    std::vector<double> v({5., -3., 10.});
    IMP::Vec<3,double> vec(v);

    // Check.
    EXPECT_DOUBLE_EQ(5., vec.At(0));
    EXPECT_DOUBLE_EQ(-3., vec.At(1));
    EXPECT_DOUBLE_EQ(10., vec.At(2));

}


// Check constructor from dynamic memory array. More extensive test in the corresponding Set() function.
TEST_F(VecTest, InitFromDynMem) {

    // Set to random values.
    double *mem = new double[3];
    mem[0] = 5.; mem[1] = -3.; mem[2] = 10.;
    IMP::Vec<3,double> vec(mem, 3);

    // Check.
    EXPECT_DOUBLE_EQ(5., vec.At(0));
    EXPECT_DOUBLE_EQ(-3., vec.At(1));
    EXPECT_DOUBLE_EQ(10., vec.At(2));

}


// Check setting vector data from initializer list.
TEST_F(VecTest, SetFromInitializerList) {

    // Set to random values.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Check.
    EXPECT_EQ(150, this->vec_u1_.At(0));
    EXPECT_EQ(-10, this->vec_i2_.At(0));
    EXPECT_EQ(35, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.12f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(-0.001f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(15.f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(12., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(-3., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(3.14, this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(-152.269, this->vec_d4_.At(3));

    // Check Set out of range throw.
    EXPECT_THROW(this->vec_u1_.Set({10,20}), std::invalid_argument);

}


// Check setting vector data from copying std::vector.
TEST_F(VecTest, SetFromCopyStdVector) {

    // Initialize std::vectors.
    std::vector<unsigned> vector_u1({150});
    std::vector<int> vector_i2({-10, 35});
    std::vector<float> vector_f3({0.12f, -0.001f, 15.f});
    std::vector<double> vector_d4({12., -3., 3.14, -152.269});

    // Set values from std::vectors.
    this->vec_u1_.Set(vector_u1);
    this->vec_i2_.Set(vector_i2);
    this->vec_f3_.Set(vector_f3);
    this->vec_d4_.Set(vector_d4);

    // Check.
    EXPECT_EQ(150, this->vec_u1_.At(0));
    EXPECT_EQ(-10, this->vec_i2_.At(0));
    EXPECT_EQ(35, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.12f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(-0.001f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(15.f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(12., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(-3., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(3.14, this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(-152.269, this->vec_d4_.At(3));

    // Check Set out of range throw.
    std::vector<unsigned> bad_vector_u2({10,20});
    EXPECT_THROW(this->vec_u1_.Set(bad_vector_u2), std::invalid_argument);

}


// Check setting vector data from moving std::vector.
TEST_F(VecTest, SetFromMoveStdVector) {

    // Initialize std::vectors.
    std::vector<unsigned> vector_u1({150});
    std::vector<int> vector_i2({-10, 35});
    std::vector<float> vector_f3({0.12f, -0.001f, 15.f});
    std::vector<double> vector_d4({12., -3., 3.14, -152.269});

    // Set values from std::vectors.
    this->vec_u1_.Set(std::move(vector_u1));
    this->vec_i2_.Set(std::move(vector_i2));
    this->vec_f3_.Set(std::move(vector_f3));
    this->vec_d4_.Set(std::move(vector_d4));

    // Check.
    EXPECT_EQ(150, this->vec_u1_.At(0));
    EXPECT_EQ(-10, this->vec_i2_.At(0));
    EXPECT_EQ(35, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.12f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(-0.001f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(15.f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(12., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(-3., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(3.14, this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(-152.269, this->vec_d4_.At(3));

    // Check Set out of range throw.
    std::vector<unsigned> bad_vector_u2({10,20});
    EXPECT_THROW(this->vec_u1_.Set(bad_vector_u2), std::invalid_argument);

}


// Check setting vector data from dynamic memory array.
TEST_F(VecTest, SetFromDynamicArray) {

    // Dynamic memory arrays with random values.
    unsigned *dyn_array_u1 = new unsigned[1];
    dyn_array_u1[0] = 150;

    int *dyn_array_i2 = new int[2];
    dyn_array_i2[0] = -10; dyn_array_i2[1] = 35;

    float *dyn_array_f3 = new float[3];
    dyn_array_f3[0] = 0.12f; dyn_array_f3[1] = -0.001f; dyn_array_f3[2] = 15.f;

    double *dyn_array_d4 = new double[4];
    dyn_array_d4[0] = 12.; dyn_array_d4[1] = -3.; dyn_array_d4[2] = 3.14; dyn_array_d4[3] = -152.269;

    // Set vector data.
    this->vec_u1_.Set(dyn_array_u1, 1);
    this->vec_i2_.Set(dyn_array_i2, 2);
    this->vec_f3_.Set(dyn_array_f3, 3);
    this->vec_d4_.Set(dyn_array_d4, 4);

    // Check.
    EXPECT_EQ(150, this->vec_u1_.At(0));
    EXPECT_EQ(-10, this->vec_i2_.At(0));
    EXPECT_EQ(35, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.12f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(-0.001f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(15.f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(12., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(-3., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(3.14, this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(-152.269, this->vec_d4_.At(3));

    // Check Set out of range throw.
    unsigned *bad_dyn_array_u2 = new unsigned[2];
    bad_dyn_array_u2[0] = 10; bad_dyn_array_u2[1] = 20;

    EXPECT_THROW(this->vec_u1_.Set(bad_dyn_array_u2, 2), std::invalid_argument);

    // Clear memory.
    delete [] dyn_array_u1;
    delete [] dyn_array_i2;
    delete [] dyn_array_f3;
    delete [] dyn_array_d4;
    delete [] bad_dyn_array_u2;

}


// Check component wise multiplication with other vector.
TEST_F(VecTest, CwiseMul) {

    // Set to random values.
    this->vec_u1_.Set({2});
    this->vec_i2_.Set({-2, 2});
    this->vec_f3_.Set({2.5f, -3.2f, 15.f});
    this->vec_d4_.Set({7., -1., 50., -2.5});

    // Create other vectors.
    IMP::Vec<1, unsigned> other_vec_u1({3});
    IMP::Vec<2, int> other_vec_i2({5, 4});
    IMP::Vec<3, float> other_vec_f3({2.f, 3.f, 1.f});
    IMP::Vec<4, double> other_vec_d4({2., -3., 5., 1.});

    // Component-wise multiplication results.
    IMP::Vec<1, unsigned> res_u1 = this->vec_u1_.CwiseMul(other_vec_u1);
    IMP::Vec<2, int> res_i2 = this->vec_i2_.CwiseMul(other_vec_i2);
    IMP::Vec<3, float> res_f3 = this->vec_f3_.CwiseMul(other_vec_f3);
    IMP::Vec<4, double> res_d4 = this->vec_d4_.CwiseMul(other_vec_d4);

    // Check.
    EXPECT_EQ(6, res_u1.At(0));
    EXPECT_EQ(-10, res_i2.At(0));
    EXPECT_EQ(8, res_i2.At(1));
    EXPECT_FLOAT_EQ(5.f, res_f3.At(0));
    EXPECT_FLOAT_EQ(-9.6f, res_f3.At(1));
    EXPECT_FLOAT_EQ(15.f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(14., res_d4.At(0));
    EXPECT_DOUBLE_EQ(3., res_d4.At(1));
    EXPECT_DOUBLE_EQ(250., res_d4.At(2));
    EXPECT_DOUBLE_EQ(-2.5, res_d4.At(3));

}


// Check component wise division by other vector.
TEST_F(VecTest, CwiseDiv) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-10, 8});
    this->vec_f3_.Set({5.f, -9.6f, 15.f});
    this->vec_d4_.Set({14., 3., 250., -2.5});

    // Create other vectors.
    IMP::Vec<1, unsigned> other_vec_u1({3});
    IMP::Vec<2, int> other_vec_i2({5, 4});
    IMP::Vec<3, float> other_vec_f3({2.f, 3.f, 1.f});
    IMP::Vec<4, double> other_vec_d4({2., -3., 5., 1.});

    // Component-wise multiplication results.
    IMP::Vec<1, unsigned> res_u1 = this->vec_u1_.CwiseDiv(other_vec_u1);
    IMP::Vec<2, int> res_i2 = this->vec_i2_.CwiseDiv(other_vec_i2);
    IMP::Vec<3, float> res_f3 = this->vec_f3_.CwiseDiv(other_vec_f3);
    IMP::Vec<4, double> res_d4 = this->vec_d4_.CwiseDiv(other_vec_d4);

    // Check.
    EXPECT_EQ(2, res_u1.At(0));
    EXPECT_EQ(-2, res_i2.At(0));
    EXPECT_EQ(2, res_i2.At(1));
    EXPECT_FLOAT_EQ(2.5f, res_f3.At(0));
    EXPECT_FLOAT_EQ(-3.2f, res_f3.At(1));
    EXPECT_FLOAT_EQ(15.f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(7., res_d4.At(0));
    EXPECT_DOUBLE_EQ(-1., res_d4.At(1));
    EXPECT_DOUBLE_EQ(50., res_d4.At(2));
    EXPECT_DOUBLE_EQ(-2.5, res_d4.At(3));

}


// Check component wise absolute value.
TEST_F(VecTest, CwiseAbs) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-10, 8});
    this->vec_f3_.Set({5.f, -9.6f, 15.f});
    this->vec_d4_.Set({14., -3., -250., -2.5});

    // Component-wise absolute value results.
    IMP::Vec<1, unsigned> res_u1 = this->vec_u1_.CwiseAbs();
    IMP::Vec<2, int> res_i2 = this->vec_i2_.CwiseAbs();
    IMP::Vec<3, float> res_f3 = this->vec_f3_.CwiseAbs();
    IMP::Vec<4, double> res_d4 = this->vec_d4_.CwiseAbs();

    // Check.
    EXPECT_EQ(6, res_u1.At(0));
    EXPECT_EQ(10, res_i2.At(0));
    EXPECT_EQ(8, res_i2.At(1));
    EXPECT_FLOAT_EQ(5.f, res_f3.At(0));
    EXPECT_FLOAT_EQ(9.6f, res_f3.At(1));
    EXPECT_FLOAT_EQ(15.f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(14., res_d4.At(0));
    EXPECT_DOUBLE_EQ(3., res_d4.At(1));
    EXPECT_DOUBLE_EQ(250., res_d4.At(2));
    EXPECT_DOUBLE_EQ(2.5, res_d4.At(3));

}


// Check squared length computation.
TEST_F(VecTest, Length2) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-2, 1});
    this->vec_f3_.Set({5.f, 2.f, 3.2f});
    this->vec_d4_.Set({3., -1., -2.5, 2.5});

    // Component-wise absolute value results.
    unsigned sq_length_u1 = this->vec_u1_.Length2();
    int sq_length_i2 = this->vec_i2_.Length2();
    float sq_length_f3 = this->vec_f3_.Length2();
    double sq_length_d4 = this->vec_d4_.Length2();

    // Check.
    EXPECT_EQ(36, sq_length_u1);
    EXPECT_EQ(5, sq_length_i2);
    EXPECT_FLOAT_EQ(39.24f, sq_length_f3);
    EXPECT_DOUBLE_EQ(22.5, sq_length_d4);

}


// Check squared distance computation.
TEST_F(VecTest, Distance2) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-2, 1});
    this->vec_f3_.Set({5.f, 2.f, 3.2f});
    this->vec_d4_.Set({3., -1., -2.5, 2.5});

    // Set other vectors values.
    IMP::Vec<1, unsigned> other_vec_u1({6});
    IMP::Vec<2, int> other_vec_i2({-2, 1});
    IMP::Vec<3, float> other_vec_f3({5.f, 2.f, 3.2f});
    IMP::Vec<4, double> other_vec_d4({3., -1., -2.5, 2.5});


    // Component-wise absolute value results.
    unsigned sq_dist_u1 = this->vec_u1_.Distance2(other_vec_u1);
    int sq_dist_i2 = this->vec_i2_.Distance2(other_vec_i2);
    float sq_dist_f3 = this->vec_f3_.Distance2(other_vec_f3);
    double sq_dist_d4 = this->vec_d4_.Distance2(other_vec_d4);

    // Check.
    EXPECT_EQ(0, sq_dist_u1);
    EXPECT_EQ(0, sq_dist_i2);
    EXPECT_FLOAT_EQ(0.f, sq_dist_f3);
    EXPECT_DOUBLE_EQ(0., sq_dist_d4);

}


// Check maximum coeffient access.
TEST_F(VecTest, MaxCoeff) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-2, 1});
    this->vec_f3_.Set({5.f, 2.f, 3.2f});
    this->vec_d4_.Set({3., -1., -2.5, 2.5});

    // Get maximum coefficients.
    unsigned max_u1 = this->vec_u1_.MaxCoeff();
    int max_i2 = this->vec_i2_.MaxCoeff();
    float max_f3 = this->vec_f3_.MaxCoeff();
    double max_d4 = this->vec_d4_.MaxCoeff();

    // Check.
    EXPECT_EQ(6, max_u1);
    EXPECT_EQ(1, max_i2);
    EXPECT_FLOAT_EQ(5.f, max_f3);
    EXPECT_DOUBLE_EQ(3., max_d4);

}


// Check minimum coeffient access.
TEST_F(VecTest, MinCoeff) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-2, 1});
    this->vec_f3_.Set({5.f, 2.f, 3.2f});
    this->vec_d4_.Set({3., -1., -2.5, 2.5});

    // Get minimum coefficients.
    unsigned min_u1 = this->vec_u1_.MinCoeff();
    int min_i2 = this->vec_i2_.MinCoeff();
    float min_f3 = this->vec_f3_.MinCoeff();
    double min_d4 = this->vec_d4_.MinCoeff();

    // Check.
    EXPECT_EQ(6, min_u1);
    EXPECT_EQ(-2, min_i2);
    EXPECT_FLOAT_EQ(2.f, min_f3);
    EXPECT_DOUBLE_EQ(-2.5, min_d4);

}


// Check set to zero.
TEST_F(VecTest, SetZero) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-2, 1});
    this->vec_f3_.Set({5.f, 2.f, 3.2f});
    this->vec_d4_.Set({3., -1., -2.5, 2.5});

    // Check.
    EXPECT_EQ(6, this->vec_u1_.At(0));
    EXPECT_EQ(-2, this->vec_i2_.At(0));
    EXPECT_EQ(1, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(5.f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(2.f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(3.2f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(3., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(-1., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(-2.5, this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(2.5, this->vec_d4_.At(3));

    // Set vectors to zero.
    this->vec_u1_.SetZero();
    this->vec_i2_.SetZero();
    this->vec_f3_.SetZero();
    this->vec_d4_.SetZero();

    // Check after setting zero.
    EXPECT_EQ(0, this->vec_u1_.At(0));
    EXPECT_EQ(0, this->vec_i2_.At(0));
    EXPECT_EQ(0, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(0.f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(0.f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(0., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(0., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(0., this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(0., this->vec_d4_.At(3));

}


// Check is zero.
TEST_F(VecTest, IsZero) {

    // Check.
    EXPECT_TRUE(this->vec_u1_.IsZero());
    EXPECT_TRUE(this->vec_i2_.IsZero());
    EXPECT_TRUE(this->vec_f3_.IsZero());
    EXPECT_TRUE(this->vec_d4_.IsZero());

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-2, 0});
    this->vec_f3_.Set({0.f, 2.f, 0.f});
    this->vec_d4_.Set({3., -1., -2.5, 2.5});


    // Check after setting values.
    EXPECT_FALSE(this->vec_u1_.IsZero());
    EXPECT_FALSE(this->vec_i2_.IsZero());
    EXPECT_FALSE(this->vec_f3_.IsZero());
    EXPECT_FALSE(this->vec_d4_.IsZero());

}


// Check range loop.
TEST_F(VecTest, RangeLoop) {

    // Set to random values.
    this->vec_u1_.Set({6});
    this->vec_i2_.Set({-2, 0});
    this->vec_f3_.Set({0.f, 2.f, 0.f});
    this->vec_d4_.Set({3., -1., -2.5, 2.5});

    // Create other vectors with same values.
    IMP::Vec<2, int> other_vec_i2; other_vec_i2 = this->vec_i2_;
    IMP::Vec<3, float> other_vec_f3; other_vec_f3 = this->vec_f3_;
    IMP::Vec<4, double> other_vec_d4; other_vec_d4 = this->vec_d4_;

    //Check loops.
    for (auto val_u : this->vec_u1_) {
        EXPECT_EQ(6, val_u);
    }
    for (auto &val_i : this->vec_i2_) {
        auto id = &val_i - &this->vec_i2_[0];
        EXPECT_EQ(val_i, other_vec_i2.At(id));
    }
    int f_id = 0;
    for (const auto val_f : this->vec_f3_) {
        EXPECT_FLOAT_EQ(val_f, other_vec_f3.At(f_id));
        f_id++;
    }
    f_id = 0;
    for (const auto &val_d : this->vec_d4_) {
        auto id = &val_d - &this->vec_d4_[0];
        EXPECT_DOUBLE_EQ(val_d, other_vec_d4.At(id));
    }

}


// Check accessing and modifying vector data with square brackets operator.
TEST_F(VecTest, SqBracketsAccess) {

    // Set to random values.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Check.
    EXPECT_EQ(150, this->vec_u1_[0]);
    EXPECT_EQ(-10, this->vec_i2_[0]);
    EXPECT_EQ(35, this->vec_i2_[1]);
    EXPECT_FLOAT_EQ(0.12f, this->vec_f3_[0]);
    EXPECT_FLOAT_EQ(-0.001f, this->vec_f3_[1]);
    EXPECT_FLOAT_EQ(15.f, this->vec_f3_[2]);
    EXPECT_DOUBLE_EQ(12., this->vec_d4_[0]);
    EXPECT_DOUBLE_EQ(-3., this->vec_d4_[1]);
    EXPECT_DOUBLE_EQ(3.14, this->vec_d4_[2]);
    EXPECT_DOUBLE_EQ(-152.269, this->vec_d4_[3]);

    // Modify values.
    this->vec_u1_[0] = 140;
    this->vec_i2_[0] = -5; this->vec_i2_[1] = 40;
    this->vec_f3_[0] = 0.f; this->vec_f3_[1] = -5.f; this->vec_f3_[2] = 10.f;
    this->vec_d4_[0] = 15.; this->vec_d4_[1] = 20.;
    this->vec_d4_[2] = -10.; this->vec_d4_[3] = -5.5;

    // Check.
    EXPECT_EQ(140, this->vec_u1_[0]);
    EXPECT_EQ(-5, this->vec_i2_[0]);
    EXPECT_EQ(40, this->vec_i2_[1]);
    EXPECT_FLOAT_EQ(0.f, this->vec_f3_[0]);
    EXPECT_FLOAT_EQ(-5.f, this->vec_f3_[1]);
    EXPECT_FLOAT_EQ(10.f, this->vec_f3_[2]);
    EXPECT_DOUBLE_EQ(15., this->vec_d4_[0]);
    EXPECT_DOUBLE_EQ(20., this->vec_d4_[1]);
    EXPECT_DOUBLE_EQ(-10., this->vec_d4_[2]);
    EXPECT_DOUBLE_EQ(-5.5, this->vec_d4_[3]);

    // Check that doesn't throw.
    EXPECT_NO_THROW(this->vec_u1_[5]);
    EXPECT_NO_THROW(this->vec_d4_[10]);

}


// Check accessing and modifying vector data with At() function.
TEST_F(VecTest, AtAccess) {

    // Set to random values.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Check.
    EXPECT_EQ(150, this->vec_u1_.At(0));
    EXPECT_EQ(-10, this->vec_i2_.At(0));
    EXPECT_EQ(35, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.12f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(-0.001f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(15.f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(12., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(-3., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(3.14, this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(-152.269, this->vec_d4_.At(3));

    // Modify values.
    this->vec_u1_.At(0) = 140;
    this->vec_i2_.At(0) = -5; this->vec_i2_.At(1) = 40;
    this->vec_f3_.At(0) = 0.f; this->vec_f3_.At(1) = -5.f; this->vec_f3_.At(2) = 10.f;
    this->vec_d4_.At(0) = 15.; this->vec_d4_.At(1) = 20.;
    this->vec_d4_.At(2) = -10.; this->vec_d4_.At(3) = -5.5;

    // Check.
    EXPECT_EQ(140, this->vec_u1_.At(0));
    EXPECT_EQ(-5, this->vec_i2_.At(0));
    EXPECT_EQ(40, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(-5.f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(10.f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(15., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(20., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(-10., this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(-5.5, this->vec_d4_.At(3));

    // Check that throws.
    EXPECT_THROW(this->vec_u1_.At(5), std::out_of_range);
    EXPECT_THROW(this->vec_d4_.At(10), std::out_of_range);

}


// Check Vec equal operator.
TEST_F(VecTest, EqualOperator) {

    // Set values of vectors.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Create other vectors with same values.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({150});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-10, 35});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.12f, -0.001f, 15.f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({12., -3., 3.14, -152.269});
    IMP::Vec<4, double> bad_vec_d4; bad_vec_d4.Set({12.1, -4., 3.1, -152.268});

    // Check.
    EXPECT_TRUE(this->vec_u1_ == other_vec_u1);
    EXPECT_TRUE(this->vec_i2_ == other_vec_i2);
    EXPECT_TRUE(this->vec_f3_ == other_vec_f3);
    EXPECT_TRUE(this->vec_d4_ == other_vec_d4);
    EXPECT_FALSE(this->vec_d4_ == bad_vec_d4);

}


// Check Vec not equal operator.
TEST_F(VecTest, NotEqualOperator) {

    // Set values of vectors.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Create other vectors with same values.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({151});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-11, 34});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.1f, -0.005f, 12.f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({12.0001, -3., 3.1454, -152.269});
    IMP::Vec<4, double> bad_vec_d4; bad_vec_d4.Set({12., -3., 3.14, -152.269});

    // Check.
    EXPECT_TRUE(this->vec_u1_ != other_vec_u1);
    EXPECT_TRUE(this->vec_i2_ != other_vec_i2);
    EXPECT_TRUE(this->vec_f3_ != other_vec_f3);
    EXPECT_TRUE(this->vec_d4_ != other_vec_d4);
    EXPECT_FALSE(this->vec_d4_ != bad_vec_d4);

}


// Check Vec less than operator.
TEST_F(VecTest, LessThanOperator) {

    // Set values of vectors.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Create other vectors with same values.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({151});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-9, 36});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.13f, 0.005f, 16.f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({12.0001, -2.5, 3.1454, -152.});
    IMP::Vec<3, float> bad_vec_f3; bad_vec_f3.Set({0.12f, -0.001f, 15.f});
    IMP::Vec<4, double> bad_vec_d4; bad_vec_d4.Set({13., -4.5, -3.1454, -152.2});

    // Check.
    EXPECT_TRUE(this->vec_u1_ < other_vec_u1);
    EXPECT_TRUE(this->vec_i2_ < other_vec_i2);
    EXPECT_TRUE(this->vec_f3_ < other_vec_f3);
    EXPECT_TRUE(this->vec_d4_ < other_vec_d4);
    EXPECT_FALSE(this->vec_f3_ < bad_vec_f3);
    EXPECT_FALSE(this->vec_d4_ < bad_vec_d4);

}


// Check Vec greater than operator.
TEST_F(VecTest, GreaterThanOperator) {

    // Set values of vectors.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Create other vectors with same values.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({140});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-15, 30});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.1f, -0.005f, 10.f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({11.05, -5., 2.5, -152.2690000001});
    IMP::Vec<3, float> bad_vec_f3; bad_vec_f3.Set({0.12f, -0.001f, 15.f});
    IMP::Vec<4, double> bad_vec_d4; bad_vec_d4.Set({13., -4.5, -3.1454, -152.2});

    // Check.
    EXPECT_TRUE(this->vec_u1_ > other_vec_u1);
    EXPECT_TRUE(this->vec_i2_ > other_vec_i2);
    EXPECT_TRUE(this->vec_f3_ > other_vec_f3);
    EXPECT_TRUE(this->vec_d4_ > other_vec_d4);
    EXPECT_FALSE(this->vec_f3_ > bad_vec_f3);
    EXPECT_FALSE(this->vec_d4_ > bad_vec_d4);

}


// Check Vec less or equal operator.
TEST_F(VecTest, LessOrEqualOperator) {

    // Set values of vectors.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Create other vectors with same values.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({151});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-10, 35});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.13f, 0.005f, 16.f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({12., -3., 3.14, -152.269});
    IMP::Vec<4, double> bad_vec_d4; bad_vec_d4.Set({12., -3., 3.14, -154.});

    // Check.
    EXPECT_TRUE(this->vec_u1_ <= other_vec_u1);
    EXPECT_TRUE(this->vec_i2_ <= other_vec_i2);
    EXPECT_TRUE(this->vec_f3_ <= other_vec_f3);
    EXPECT_TRUE(this->vec_d4_ <= other_vec_d4);
    EXPECT_FALSE(this->vec_d4_ <= bad_vec_d4);

}


// Check Vec greater or equal operator.
TEST_F(VecTest, GreaterOrEqualOperator) {

    // Set values of vectors.
    this->vec_u1_.Set({150});
    this->vec_i2_.Set({-10, 35});
    this->vec_f3_.Set({0.12f, -0.001f, 15.f});
    this->vec_d4_.Set({12., -3., 3.14, -152.269});

    // Create other vectors with same values.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({149});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-10, 35});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.11f, -0.0015f, 14.9f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({12., -3., 3.14, -152.269});
    IMP::Vec<4, double> bad_vec_d4; bad_vec_d4.Set({12., -3., 3.14, -152.268});

    // Check.
    EXPECT_TRUE(this->vec_u1_ >= other_vec_u1);
    EXPECT_TRUE(this->vec_i2_ >= other_vec_i2);
    EXPECT_TRUE(this->vec_f3_ >= other_vec_f3);
    EXPECT_TRUE(this->vec_d4_ >= other_vec_d4);
    EXPECT_FALSE(this->vec_d4_ >= bad_vec_d4);

}


// Check Vec copy assign operator.
TEST_F(VecTest, CopyAssignOperator) {

    // Create other vectors.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({149});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-10, 35});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.11f, -0.0015f, 14.9f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({12., -3., 3.14, -152.269});

    // Assign values from other vectors.
    this->vec_u1_ = other_vec_u1;
    this->vec_i2_ = other_vec_i2;
    this->vec_f3_ = other_vec_f3;
    this->vec_d4_ = other_vec_d4;

    // Check.
    EXPECT_TRUE(this->vec_u1_ == other_vec_u1);
    EXPECT_TRUE(this->vec_i2_ == other_vec_i2);
    EXPECT_TRUE(this->vec_f3_ == other_vec_f3);
    EXPECT_TRUE(this->vec_d4_ == other_vec_d4);

}


// Check Vec move assign operator.
TEST_F(VecTest, MoveAssignOperator) {

    // Create other vectors.
    IMP::Vec<1, unsigned> other_vec_u1; other_vec_u1.Set({149});
    IMP::Vec<2, int> other_vec_i2; other_vec_i2.Set({-10, 35});
    IMP::Vec<3, float> other_vec_f3; other_vec_f3.Set({0.11f, -0.0015f, 14.9f});
    IMP::Vec<4, double> other_vec_d4; other_vec_d4.Set({12., -3., 3.14, -152.269});

    // Assign values from other vectors.
    this->vec_u1_ = std::move(other_vec_u1);
    this->vec_i2_ = std::move(other_vec_i2);
    this->vec_f3_ = std::move(other_vec_f3);
    this->vec_d4_ = std::move(other_vec_d4);

    // Check.
    EXPECT_EQ(149, this->vec_u1_.At(0));
    EXPECT_EQ(-10, this->vec_i2_.At(0));
    EXPECT_EQ(35, this->vec_i2_.At(1));
    EXPECT_FLOAT_EQ(0.11f, this->vec_f3_.At(0));
    EXPECT_FLOAT_EQ(-0.0015f, this->vec_f3_.At(1));
    EXPECT_FLOAT_EQ(14.9f, this->vec_f3_.At(2));
    EXPECT_DOUBLE_EQ(12., this->vec_d4_.At(0));
    EXPECT_DOUBLE_EQ(-3., this->vec_d4_.At(1));
    EXPECT_DOUBLE_EQ(3.14, this->vec_d4_.At(2));
    EXPECT_DOUBLE_EQ(-152.269, this->vec_d4_.At(3));

    EXPECT_EQ(0, std::distance(other_vec_u1.begin(), other_vec_u1.end()));
    EXPECT_EQ(0, std::distance(other_vec_i2.begin(), other_vec_i2.end()));
    EXPECT_EQ(0, std::distance(other_vec_f3.begin(), other_vec_f3.end()));
    EXPECT_EQ(0, std::distance(other_vec_d4.begin(), other_vec_d4.end()));

}


// Check Vec addition assign operator.
TEST_F(VecTest, AddAssignOperator) {

    // Create other vectors.
    IMP::Vec<1, unsigned> vec1_u1, vec2_u1;
    vec1_u1.Set({1}); vec2_u1.Set({2});

    IMP::Vec<2, int> vec1_i2, vec2_i2;
    vec1_i2.Set({-1, 1}); vec2_i2.Set({2, 3});

    IMP::Vec<3, float> vec1_f3, vec2_f3;
    vec1_f3.Set({-1.5f, 1.f, 1.5f}); vec2_f3.Set({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec1_d4, vec2_d4;
    vec1_d4.Set({0., 2., 5., -10.}); vec2_d4.Set({-3., 0., 0., 10.});

    // Add vector in place.
    vec1_u1 += vec2_u1;
    vec1_i2 += vec2_i2;
    vec1_f3 += vec2_f3;
    vec1_d4 += vec2_d4;

    // Check.
    EXPECT_EQ(3, vec1_u1.At(0));
    EXPECT_EQ(1, vec1_i2.At(0));
    EXPECT_EQ(4, vec1_i2.At(1));
    EXPECT_FLOAT_EQ(0.5f, vec1_f3.At(0));
    EXPECT_FLOAT_EQ(4.f, vec1_f3.At(1));
    EXPECT_FLOAT_EQ(-4.5f, vec1_f3.At(2));
    EXPECT_DOUBLE_EQ(-3., vec1_d4.At(0));
    EXPECT_DOUBLE_EQ(2., vec1_d4.At(1));
    EXPECT_DOUBLE_EQ(5., vec1_d4.At(2));
    EXPECT_DOUBLE_EQ(0., vec1_d4.At(3));

}


// Check Vec subtraction assign operator.
TEST_F(VecTest, SubAssignOperator) {

    // Create other vectors.
    IMP::Vec<1, unsigned> vec1_u1, vec2_u1;
    vec1_u1.Set({1}); vec2_u1.Set({1});

    IMP::Vec<2, int> vec1_i2, vec2_i2;
    vec1_i2.Set({-1, 1}); vec2_i2.Set({2, 3});

    IMP::Vec<3, float> vec1_f3, vec2_f3;
    vec1_f3.Set({-1.5f, 1.f, 1.5f}); vec2_f3.Set({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec1_d4, vec2_d4;
    vec1_d4.Set({0., 2., 5., -10.}); vec2_d4.Set({-3., 0., 0., 10.});

    // Subtract vector in place.
    vec1_u1 -= vec2_u1;
    vec1_i2 -= vec2_i2;
    vec1_f3 -= vec2_f3;
    vec1_d4 -= vec2_d4;

    // Check.
    EXPECT_EQ(0, vec1_u1.At(0));
    EXPECT_EQ(-3, vec1_i2.At(0));
    EXPECT_EQ(-2, vec1_i2.At(1));
    EXPECT_FLOAT_EQ(-3.5f, vec1_f3.At(0));
    EXPECT_FLOAT_EQ(-2.f, vec1_f3.At(1));
    EXPECT_FLOAT_EQ(7.5f, vec1_f3.At(2));
    EXPECT_DOUBLE_EQ(3., vec1_d4.At(0));
    EXPECT_DOUBLE_EQ(2., vec1_d4.At(1));
    EXPECT_DOUBLE_EQ(5., vec1_d4.At(2));
    EXPECT_DOUBLE_EQ(-20., vec1_d4.At(3));

}


// Check Vec scalar multiplication assign operator.
TEST_F(VecTest, MulAssignOperator) {

    // Create other vector.
    IMP::Vec<1, unsigned> vec_u1({1});

    IMP::Vec<2, int> vec_i2({2, 3});

    IMP::Vec<3, float> vec_f3({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec_d4({-3., 5., 0., 10.});

    // Multiply by scalar in place.
    vec_u1 *= 2;
    vec_i2 *= -2;
    vec_f3 *= 0.5f;
    vec_d4 *= -0.5;

    // Check.
    EXPECT_EQ(2, vec_u1.At(0));
    EXPECT_EQ(-4, vec_i2.At(0));
    EXPECT_EQ(-6, vec_i2.At(1));
    EXPECT_FLOAT_EQ(1.f, vec_f3.At(0));
    EXPECT_FLOAT_EQ(1.5f, vec_f3.At(1));
    EXPECT_FLOAT_EQ(-3.f, vec_f3.At(2));
    EXPECT_DOUBLE_EQ(1.5, vec_d4.At(0));
    EXPECT_DOUBLE_EQ(-2.5, vec_d4.At(1));
    EXPECT_DOUBLE_EQ(0., vec_d4.At(2));
    EXPECT_DOUBLE_EQ(-5., vec_d4.At(3));

}


// Check Vec scalar division assign operator.
TEST_F(VecTest, DivAssignOperator) {

    // Create other vector.
    IMP::Vec<1, unsigned> vec_u1({4});

    IMP::Vec<2, int> vec_i2({2, 3});

    IMP::Vec<3, float> vec_f3({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec_d4({-3., 5., 0., 10.});

    // Divide by scalar in place.
    vec_u1 /= 2;
    vec_i2 /= -2;
    vec_f3 /= 0.5f;
    vec_d4 /= -0.5;

    // Check.
    EXPECT_EQ(2, vec_u1.At(0));
    EXPECT_EQ(-1, vec_i2.At(0));
    EXPECT_EQ(-1, vec_i2.At(1));
    EXPECT_FLOAT_EQ(4.f, vec_f3.At(0));
    EXPECT_FLOAT_EQ(6.f, vec_f3.At(1));
    EXPECT_FLOAT_EQ(-12.f, vec_f3.At(2));
    EXPECT_DOUBLE_EQ(6., vec_d4.At(0));
    EXPECT_DOUBLE_EQ(-10., vec_d4.At(1));
    EXPECT_DOUBLE_EQ(0., vec_d4.At(2));
    EXPECT_DOUBLE_EQ(-20., vec_d4.At(3));

    // Check exception behavior
    EXPECT_THROW(vec_u1 /= 0, std::overflow_error);
    EXPECT_THROW(vec_d4 /= 0., std::overflow_error);

}


// Check Vec addition operator.
TEST_F(VecTest, AddOperator) {

    // Create other vectors.
    IMP::Vec<1, unsigned> vec1_u1, vec2_u1;
    vec1_u1.Set({1}); vec2_u1.Set({2});

    IMP::Vec<2, int> vec1_i2, vec2_i2;
    vec1_i2.Set({-1, 1}); vec2_i2.Set({2, 3});

    IMP::Vec<3, float> vec1_f3, vec2_f3;
    vec1_f3.Set({-1.5f, 1.f, 1.5f}); vec2_f3.Set({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec1_d4, vec2_d4;
    vec1_d4.Set({0., 2., 5., -10.}); vec2_d4.Set({-3., 0., 0., 10.});

    // Add vectors.
    IMP::Vec<1, unsigned> res_u1 = vec1_u1 + vec2_u1;
    IMP::Vec<2, int> res_i2 = vec1_i2 + vec2_i2;
    IMP::Vec<3, float> res_f3 = vec1_f3 + vec2_f3;
    IMP::Vec<4, double> res_d4 = vec1_d4 + vec2_d4;

    // Check.
    EXPECT_EQ(3, res_u1.At(0));
    EXPECT_EQ(1, res_i2.At(0));
    EXPECT_EQ(4, res_i2.At(1));
    EXPECT_FLOAT_EQ(0.5f, res_f3.At(0));
    EXPECT_FLOAT_EQ(4.f, res_f3.At(1));
    EXPECT_FLOAT_EQ(-4.5f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(-3., res_d4.At(0));
    EXPECT_DOUBLE_EQ(2., res_d4.At(1));
    EXPECT_DOUBLE_EQ(5., res_d4.At(2));
    EXPECT_DOUBLE_EQ(0., res_d4.At(3));

}


// Check Vec subtraction assign operator.
TEST_F(VecTest, SubOperator) {

    // Create other vectors.
    IMP::Vec<1, unsigned> vec1_u1, vec2_u1;
    vec1_u1.Set({1}); vec2_u1.Set({1});

    IMP::Vec<2, int> vec1_i2, vec2_i2;
    vec1_i2.Set({-1, 1}); vec2_i2.Set({2, 3});

    IMP::Vec<3, float> vec1_f3, vec2_f3;
    vec1_f3.Set({-1.5f, 1.f, 1.5f}); vec2_f3.Set({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec1_d4, vec2_d4;
    vec1_d4.Set({0., 2., 5., -10.}); vec2_d4.Set({-3., 0., 0., 10.});

    // Subtract vectors.
    IMP::Vec<1, unsigned> res_u1 = vec1_u1 - vec2_u1;
    IMP::Vec<2, int> res_i2 = vec1_i2 - vec2_i2;
    IMP::Vec<3, float> res_f3 = vec1_f3 - vec2_f3;
    IMP::Vec<4, double> res_d4 = vec1_d4 - vec2_d4;

    // Check.
    EXPECT_EQ(0, res_u1.At(0));
    EXPECT_EQ(-3, res_i2.At(0));
    EXPECT_EQ(-2, res_i2.At(1));
    EXPECT_FLOAT_EQ(-3.5f, res_f3.At(0));
    EXPECT_FLOAT_EQ(-2.f, res_f3.At(1));
    EXPECT_FLOAT_EQ(7.5f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(3., res_d4.At(0));
    EXPECT_DOUBLE_EQ(2., res_d4.At(1));
    EXPECT_DOUBLE_EQ(5., res_d4.At(2));
    EXPECT_DOUBLE_EQ(-20., res_d4.At(3));

}


// Check Vec scalar multiplication operators, both vec*scalar and scalar*vec.
TEST_F(VecTest, MulOperator) {

    // Create other vector.
    IMP::Vec<1, unsigned> vec_u1({1});

    IMP::Vec<2, int> vec_i2({2, 3});

    IMP::Vec<3, float> vec_f3({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec_d4({-3., 5., 0., 10.});

    // Multiply by vector by scalar.
    IMP::Vec<1, unsigned> res_u1 = vec_u1 * 2u;
    IMP::Vec<2, int> res_i2 = vec_i2 * (-2);
    IMP::Vec<3, float> res_f3 = vec_f3 * 0.5f;
    IMP::Vec<4, double> res_d4 = vec_d4 * (-0.5);

    // Check.
    EXPECT_EQ(2, res_u1.At(0));
    EXPECT_EQ(-4, res_i2.At(0));
    EXPECT_EQ(-6, res_i2.At(1));
    EXPECT_FLOAT_EQ(1.f, res_f3.At(0));
    EXPECT_FLOAT_EQ(1.5f, res_f3.At(1));
    EXPECT_FLOAT_EQ(-3.f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(1.5, res_d4.At(0));
    EXPECT_DOUBLE_EQ(-2.5, res_d4.At(1));
    EXPECT_DOUBLE_EQ(0., res_d4.At(2));
    EXPECT_DOUBLE_EQ(-5., res_d4.At(3));

    // Multiply by scalar with vector now.
    res_u1 = 2u * vec_u1;
    res_i2 = -2 * vec_i2;
    res_f3 = 0.5f * vec_f3;
    res_d4 = -0.5 * vec_d4;

    // Check.
    EXPECT_EQ(2, res_u1.At(0));
    EXPECT_EQ(-4, res_i2.At(0));
    EXPECT_EQ(-6, res_i2.At(1));
    EXPECT_FLOAT_EQ(1.f, res_f3.At(0));
    EXPECT_FLOAT_EQ(1.5f, res_f3.At(1));
    EXPECT_FLOAT_EQ(-3.f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(1.5, res_d4.At(0));
    EXPECT_DOUBLE_EQ(-2.5, res_d4.At(1));
    EXPECT_DOUBLE_EQ(0., res_d4.At(2));
    EXPECT_DOUBLE_EQ(-5., res_d4.At(3));

}


// Check Vec scalar division operator.
TEST_F(VecTest, DivOperator) {

    // Create other vector.
    IMP::Vec<1, unsigned> vec_u1({4});

    IMP::Vec<2, int> vec_i2({2, 3});

    IMP::Vec<3, float> vec_f3({2.f, 3.f, -6.f});

    IMP::Vec<4, double> vec_d4({-3., 5., 0., 10.});

    // Divide by scalar.
    IMP::Vec<1, unsigned> res_u1 = vec_u1 / 2u;
    IMP::Vec<2, int> res_i2 = vec_i2 / (-2);
    IMP::Vec<3, float> res_f3 = vec_f3 / 0.5f;
    IMP::Vec<4, double> res_d4 = vec_d4 / (-0.5);

    // Check.
    EXPECT_EQ(2, res_u1.At(0));
    EXPECT_EQ(-1, res_i2.At(0));
    EXPECT_EQ(-1, res_i2.At(1));
    EXPECT_FLOAT_EQ(4.f, res_f3.At(0));
    EXPECT_FLOAT_EQ(6.f, res_f3.At(1));
    EXPECT_FLOAT_EQ(-12.f, res_f3.At(2));
    EXPECT_DOUBLE_EQ(6., res_d4.At(0));
    EXPECT_DOUBLE_EQ(-10., res_d4.At(1));
    EXPECT_DOUBLE_EQ(0., res_d4.At(2));
    EXPECT_DOUBLE_EQ(-20., res_d4.At(3));

    // Check exception behavior
    EXPECT_THROW(vec_u1 / 0u, std::overflow_error);
    EXPECT_THROW(vec_d4 / 0., std::overflow_error);

}


} // End of namespace IMP_test.
