#include "half_facet_test.hpp"

namespace IMP_test {


// Check default constructor with data initialization.
TEST_F(HalfFacetTest, InitCtor) {

    // Construct new half-facet
    IMP::HalfFacet hf(5,10,3);

    // Check if data are initialized as expected.
    EXPECT_EQ(5, hf.VertexId());
    EXPECT_EQ(10, hf.CellId());
    EXPECT_EQ(3, hf.FacetId());

}


// Check default constructor with no data initialization.
TEST_F(HalfFacetTest, DefaultCtor) {

    // Construct new half-facet
    IMP::HalfFacet hf;

    // Check if data are not initialized as expected.
    EXPECT_EQ(-1, hf.VertexId());
    EXPECT_EQ(-1, hf.CellId());
    EXPECT_EQ(-1, hf.FacetId());

}


// Check Set function.
TEST_F(HalfFacetTest, Set) {

    // Construct new half-facet
    IMP::HalfFacet hf;

    // Check if data are set correctly.
    hf.Set(3, 12, 5);
    EXPECT_EQ(3, hf.VertexId());
    EXPECT_EQ(12, hf.CellId());
    EXPECT_EQ(5, hf.FacetId());

    // Check if data are set correctly.
    hf.Set(0, 0, 0);
    EXPECT_EQ(0, hf.VertexId());
    EXPECT_EQ(0, hf.CellId());
    EXPECT_EQ(0, hf.FacetId());

    // Check if throws when negative data is given.
    EXPECT_THROW(hf.Set(-12, 0, 0), std::invalid_argument);
    EXPECT_THROW(hf.Set(0, -5, 0), std::invalid_argument);
    EXPECT_THROW(hf.Set(0, 0, -3), std::invalid_argument);
    EXPECT_THROW(hf.Set(-12, -5, -3), std::invalid_argument);

}


// Check Set function by copying other vertex-half-facet's data.
TEST_F(HalfFacetTest, SetByCopy) {

    // Construct new half-facets
    IMP::HalfFacet hf1, hf2;

    // Set hf1 data.
    hf1.Set(2, 35, 42);

    // Check if hf2 data are set correctly.
    hf2.Set(hf1);
    EXPECT_EQ(hf1.VertexId(), hf2.VertexId());
    EXPECT_EQ(2, hf2.VertexId());

    EXPECT_EQ(hf1.CellId(), hf2.CellId());
    EXPECT_EQ(35, hf2.CellId());

    EXPECT_EQ(hf1.FacetId(), hf2.FacetId());
    EXPECT_EQ(42, hf2.FacetId());

}


// Check IsInitialized function.
TEST_F(HalfFacetTest, IsInitialized) {

    // Construct new half-facet
    IMP::HalfFacet hf;

    // Check if it is not initialized yet.
    EXPECT_TRUE(hf.IsNull());

    // Set hf1 data.
    hf.Set(5,100,50);

    // Should be initialized now.
    EXPECT_FALSE(hf.IsNull());

}


// Check HalfFacet equal operator.
TEST_F(HalfFacetTest, EqualOperator) {

    // Construct new vertex-half-facets
    IMP::HalfFacet hf1, hf2;

    // Set vhf1 data.
    hf1.Set(4,35,42);

    // Set hf2 data.
    hf2.Set(3,15,27);

    // Chech if equal. Should not.
    EXPECT_FALSE(hf1 == hf2);

    // Reset hf2 data.
    hf2.Set(hf1);

    // Check if equal. Now they should.
    EXPECT_TRUE(hf1 == hf2);
    EXPECT_TRUE(hf2 == hf1);

}



} // End of namespace IMP_test.
