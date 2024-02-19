#include "cell_test.hpp"

namespace IMP_test {


// Check default constructor.
TEST_F(CellTest, InitCtor) {

    // Construct new ahf cell.
    IMP::Cell<3, 4> cell;

    // Check if data are initialized as expected.
    EXPECT_EQ(cell.SibHfacets().size(), 4);
    EXPECT_EQ(cell.SibHfacets().capacity(), 4);
    EXPECT_EQ(cell.Connectivity().Dim(), 4);
    EXPECT_TRUE(cell.Connectivity().IsZero());
}


// Check SetConnectivity function for IMP::Vec input.
TEST_F(CellTest, SetConnectivityByImpVec) {

    // Construct new 2D ahf cell.
    IMP::Cell<2,4> cell_2d;

    // IMP vector defining the 2d connectivity.
    IMP::Vec<4, int> conn_2d({15,45,0,795});

    // Set cell_2d connectivity.
    cell_2d.SetConnectivity(conn_2d);

    // Check the assigned connectivity.
    EXPECT_EQ(cell_2d.Connectivity()[0], 15);
    EXPECT_EQ(cell_2d.Connectivity()[1], 45);
    EXPECT_EQ(cell_2d.Connectivity()[2], 0);
    EXPECT_EQ(cell_2d.Connectivity()[3], 795);

    // Construct new 3D ahf cell.
    IMP::Cell<3,4> cell_3d;

    // IMP vector defining the 3d connectivity.
    IMP::Vec<4, int> conn_3d({15,45,0,795});

    // Set cell_3d connectivity.
    cell_3d.SetConnectivity(conn_3d);

    // Check the assigned connectivity.
    EXPECT_EQ(cell_3d.Connectivity()[0], 15);
    EXPECT_EQ(cell_3d.Connectivity()[1], 45);
    EXPECT_EQ(cell_3d.Connectivity()[2], 0);
    EXPECT_EQ(cell_3d.Connectivity()[3], 795);

}


// Check SetConnectivity function for std::initializer_list input.
TEST_F(CellTest, SetConnectivityByInitList) {

    // Construct new 2D ahf cell.
    IMP::Cell<2,4> cell_2d;

    // Set cell_2d connectivity.
    cell_2d.SetConnectivity({15,45,0,795});

    // Check the assigned connectivity.
    EXPECT_EQ(cell_2d.Connectivity()[0], 15);
    EXPECT_EQ(cell_2d.Connectivity()[1], 45);
    EXPECT_EQ(cell_2d.Connectivity()[2], 0);
    EXPECT_EQ(cell_2d.Connectivity()[3], 795);

    // Construct new 3D ahf cell.
    IMP::Cell<3,4> cell_3d;

    // Set cell_3d connectivity.
    cell_3d.SetConnectivity({15,45,0,795});

    // Check the assigned connectivity.
    EXPECT_EQ(cell_3d.Connectivity()[0], 15);
    EXPECT_EQ(cell_3d.Connectivity()[1], 45);
    EXPECT_EQ(cell_3d.Connectivity()[2], 0);
    EXPECT_EQ(cell_3d.Connectivity()[3], 795);

}


// Check SetAllSiblingHalfFacets function.
TEST_F(CellTest, SetAllSiblingHalfFacets) {

    // Construct new ahf cell.
    IMP::Cell<3,2> cell;

    // Construct two halfacets.
    IMP::HalfFacet hf1, hf2;

    // Set half facets data.
    hf1.Set(150,525,3);
    hf2.Set(180,525,2);

    // Store half facets in std::vector.
    std::vector<IMP::HalfFacet> hfacets{hf1, hf2};

    // Set half facets of the ahf cell.
    cell.SetAllSibHfacets(hfacets);

    // Check the assigned half facets.
    EXPECT_EQ(cell.SibHfacets()[0].VertexId(), 150);
    EXPECT_EQ(cell.SibHfacets()[0].CellId(), 525);
    EXPECT_EQ(cell.SibHfacets()[0].FacetId(), 3);

    EXPECT_EQ(cell.SibHfacets()[1].VertexId(), 180);
    EXPECT_EQ(cell.SibHfacets()[1].CellId(), 525);
    EXPECT_EQ(cell.SibHfacets()[1].FacetId(), 2);

    // Store extra half facet in std::vector.
    hfacets.emplace_back(IMP::HalfFacet(10,525,0));

    // Check if throws as larger number of half facets
    // than the cell's vertices was requested.
    EXPECT_THROW(cell.SetAllSibHfacets(hfacets), std::invalid_argument);

}


// Check SetSiblingHalfFacet function.
TEST_F(CellTest, SetSiblingHalfFacet) {

    // Construct new ahf cell.
    IMP::Cell<3,2> cell;

    // Construct two halfacets.
    IMP::HalfFacet hf1, hf2;

    // Set half facets data.
    hf1.Set(150,525,3);
    hf2.Set(180,525,2);

    // Set half facets of the ahf cell.
    cell.SetSibHfacet(0, hf1);
    cell.SetSibHfacet(1, hf2);

    // Check the assigned half facets.
    EXPECT_EQ(cell.SibHfacets()[0].VertexId(), 150);
    EXPECT_EQ(cell.SibHfacets()[0].CellId(), 525);
    EXPECT_EQ(cell.SibHfacets()[0].FacetId(), 3);

    EXPECT_EQ(cell.SibHfacets()[1].VertexId(), 180);
    EXPECT_EQ(cell.SibHfacets()[1].CellId(), 525);
    EXPECT_EQ(cell.SibHfacets()[1].FacetId(), 2);

    // Check if throws as expected for larger number
    // of half facets than the cells vertices.
    EXPECT_THROW(cell.SetSibHfacet(2, IMP::HalfFacet(10,525,0)), std::invalid_argument);

}


// Check Connectivity function.
TEST_F(CellTest, Connectivity) {

    // Construct new ahf cell.
    IMP::Cell<3,4> cell;

    // Set cell's connectivity.
    cell.SetConnectivity({5,10,98,0});

    // Check the assigned connectivity.
    EXPECT_EQ(cell.Connectivity()[0], 5);
    EXPECT_EQ(cell.Connectivity()[1], 10);
    EXPECT_EQ(cell.Connectivity()[2], 98);
    EXPECT_EQ(cell.Connectivity()[3], 0);

}

// Check Cell equal operator.
TEST_F(CellTest, EqualOperator) {

    // Construct new ahf cell.
    IMP::Cell<3,3> c1;
    c1.SetConnectivity({1,4,3});
    std::vector<IMP::HalfFacet> c1_hfs{IMP::HalfFacet(2,0,8),
                                       IMP::HalfFacet(1,0,5),
                                       IMP::HalfFacet(3,0,7)};
    c1.SetAllSibHfacets(c1_hfs);

    // Construct new ahf cell.
    IMP::Cell<3,3> c2;
    c2.SetConnectivity({1,4,3});
    std::vector<IMP::HalfFacet> c2_hfs{IMP::HalfFacet(2,0,8),
                                       IMP::HalfFacet(1,0,5),
                                       IMP::HalfFacet(3,0,7)};
    c2.SetAllSibHfacets(c2_hfs);

    // Construct new ahf cell.
    IMP::Cell<3,3> c3;
    c3.SetConnectivity({1,5,3});
    std::vector<IMP::HalfFacet> c3_hfs{IMP::HalfFacet(2,0,8),
                                       IMP::HalfFacet(1,0,5),
                                       IMP::HalfFacet(3,0,7)};
    c3.SetAllSibHfacets(c3_hfs);

    // Construct new ahf cell.
    IMP::Cell<3,3> c4;
    c4.SetConnectivity({1,4,3});
    std::vector<IMP::HalfFacet> c4_hfs{IMP::HalfFacet(2,0,8),
                                       IMP::HalfFacet(7,0,5),
                                       IMP::HalfFacet(3,0,7)};
    c4.SetAllSibHfacets(c4_hfs);

    // Chech equality.
    EXPECT_TRUE(c1 == c2);
    EXPECT_FALSE(c1 == c3);
    EXPECT_FALSE(c1 == c4);

}



} // End of namespace IMP_test.
