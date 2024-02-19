#include "mesh_test.hpp"

namespace IMP_test {


// Check default constructor.
TEST_F(MeshTest, Ctor) {

    // Construct new mesh.
    IMP::Mesh<3, 4> mesh;

    // Check if data are initialized as expected.
    EXPECT_EQ(mesh.NodesNum(), 0);
    EXPECT_EQ(mesh.CellsNum(), 0);
    EXPECT_EQ(mesh.NodesToHalfFacets().Size(), 0);
    EXPECT_FALSE(mesh.IsAhfBuilt());

}


// Check AppendNode.
TEST_F(MeshTest, AppendNode) {

    // Create a mesh in 2d.
    IMP::Mesh<2,3> mesh2d;

    // Create a mesh in 3d.
    IMP::Mesh<3,4> mesh3d;

    // Check that no nodes exist.
    EXPECT_EQ(0, mesh2d.NodesNum());
    EXPECT_EQ(0, mesh3d.NodesNum());

    // Append a node in the mesh2d.
    mesh2d.AppendNode({5.,2.});
    
    // Check if node was appended correctly.
    EXPECT_EQ(1, mesh2d.NodesNum());
    EXPECT_DOUBLE_EQ(5., mesh2d.Nodes()[0][0]);
    EXPECT_DOUBLE_EQ(2., mesh2d.Nodes()[0][1]);

    // Check exception when trying to append node of different dimensions.
    EXPECT_THROW(mesh2d.AppendNode({3.}), std::invalid_argument);
    EXPECT_THROW(mesh2d.AppendNode({-2.,20,45}), std::invalid_argument);

    // Append a node in the mesh3d.
    mesh3d.AppendNode({-4,2.3,100});
    EXPECT_EQ(1, mesh3d.NodesNum());
    EXPECT_DOUBLE_EQ(-4., mesh3d.Nodes()[0][0]);
    EXPECT_DOUBLE_EQ(2.3, mesh3d.Nodes()[0][1]);
    EXPECT_DOUBLE_EQ(100, mesh3d.Nodes()[0][2]);
}


// Check AppendCell.
TEST_F(MeshTest, AppendCell) {

    // Create a mesh in 2d.
    IMP::Mesh<2,3> mesh2d;

    // Create a mesh in 3d.
    IMP::Mesh<3,4> mesh3d;

    // Check initial cells number.
    EXPECT_EQ(0, mesh2d.CellsNum());
    EXPECT_EQ(0, mesh3d.CellsNum());

    // Append cells in mesh2d.
    mesh2d.AppendCell({1,5,2});
    mesh3d.AppendCell({10,4,5,0});

    // Check again cells number.
    EXPECT_EQ(1, mesh2d.CellsNum());
    EXPECT_EQ(1, mesh3d.CellsNum());

    // Check cells connectivity.
    EXPECT_EQ(1, mesh2d.Cells()[0].Connectivity()[0]);
    EXPECT_EQ(5, mesh2d.Cells()[0].Connectivity()[1]);
    EXPECT_EQ(2, mesh2d.Cells()[0].Connectivity()[2]);

    EXPECT_EQ(10, mesh3d.Cells()[0].Connectivity()[0]);
    EXPECT_EQ(4, mesh3d.Cells()[0].Connectivity()[1]);
    EXPECT_EQ(5, mesh3d.Cells()[0].Connectivity()[2]);
    EXPECT_EQ(0, mesh3d.Cells()[0].Connectivity()[3]);
}


// Check BuildAhf for 2d triangular mesh.
TEST_F(MeshTest, BuildAhfFor2dTrimesh) {

    // Create 2d mesh.
    IMP::Mesh<2,3> mesh;
    
    // Append cells in the mesh.
    mesh.AppendCell({0,1,4});
    mesh.AppendCell({1,2,4});
    mesh.AppendCell({2,3,4});
    mesh.AppendCell({3,0,4});
    
    // Create the AHF structure of the mesh.
    mesh.BuildAhf();

    // Check the sibling half facets for all the cells.
    // Cell 0 | Siblings: 4:<1,1> 4: <3,0>, -1: <-,->
    // Half facet 0:
    EXPECT_EQ(4, mesh.Cells()[0].SibHfacets()[0].VertexId());
    EXPECT_EQ(1, mesh.Cells()[0].SibHfacets()[0].CellId());
    EXPECT_EQ(1, mesh.Cells()[0].SibHfacets()[0].FacetId());
    // Half facet 1:
    EXPECT_EQ(4, mesh.Cells()[0].SibHfacets()[1].VertexId());
    EXPECT_EQ(3, mesh.Cells()[0].SibHfacets()[1].CellId());
    EXPECT_EQ(0, mesh.Cells()[0].SibHfacets()[1].FacetId());
    // Half facet 2:
    EXPECT_EQ(-1, mesh.Cells()[0].SibHfacets()[2].VertexId());
    EXPECT_EQ(-1, mesh.Cells()[0].SibHfacets()[2].CellId());
    EXPECT_EQ(-1, mesh.Cells()[0].SibHfacets()[2].FacetId());

    // Cell 1 | Siblings: 4:<2,1> 4: <0,0>, -1: <-,->
    // Half facet 0:
    EXPECT_EQ(4, mesh.Cells()[1].SibHfacets()[0].VertexId());
    EXPECT_EQ(2, mesh.Cells()[1].SibHfacets()[0].CellId());
    EXPECT_EQ(1, mesh.Cells()[1].SibHfacets()[0].FacetId());
    // Half facet 1:
    EXPECT_EQ(4, mesh.Cells()[1].SibHfacets()[1].VertexId());
    EXPECT_EQ(0, mesh.Cells()[1].SibHfacets()[1].CellId());
    EXPECT_EQ(0, mesh.Cells()[1].SibHfacets()[1].FacetId());
    // Half facet 2:
    EXPECT_EQ(-1, mesh.Cells()[1].SibHfacets()[2].VertexId());
    EXPECT_EQ(-1, mesh.Cells()[1].SibHfacets()[2].CellId());
    EXPECT_EQ(-1, mesh.Cells()[1].SibHfacets()[2].FacetId());

    // Cell 2 | Siblings: 4:<3,1> 4: <1,0>, -1: <-,->
    // Half facet 0:
    EXPECT_EQ(4, mesh.Cells()[2].SibHfacets()[0].VertexId());
    EXPECT_EQ(3, mesh.Cells()[2].SibHfacets()[0].CellId());
    EXPECT_EQ(1, mesh.Cells()[2].SibHfacets()[0].FacetId());
    // Half facet 1:
    EXPECT_EQ(4, mesh.Cells()[2].SibHfacets()[1].VertexId());
    EXPECT_EQ(1, mesh.Cells()[2].SibHfacets()[1].CellId());
    EXPECT_EQ(0, mesh.Cells()[2].SibHfacets()[1].FacetId());
    // Half facet 2:
    EXPECT_EQ(-1, mesh.Cells()[2].SibHfacets()[2].VertexId());
    EXPECT_EQ(-1, mesh.Cells()[2].SibHfacets()[2].CellId());
    EXPECT_EQ(-1, mesh.Cells()[2].SibHfacets()[2].FacetId());

    // Cell 3 | Siblings: 4:<0,1> 4: <2,0>, -1: <-,->
    // Half facet 0:
    EXPECT_EQ(4, mesh.Cells()[3].SibHfacets()[0].VertexId());
    EXPECT_EQ(0, mesh.Cells()[3].SibHfacets()[0].CellId());
    EXPECT_EQ(1, mesh.Cells()[3].SibHfacets()[0].FacetId());
    // Half facet 1:
    EXPECT_EQ(4, mesh.Cells()[3].SibHfacets()[1].VertexId());
    EXPECT_EQ(2, mesh.Cells()[3].SibHfacets()[1].CellId());
    EXPECT_EQ(0, mesh.Cells()[3].SibHfacets()[1].FacetId());
    // Half facet 2:
    EXPECT_EQ(-1, mesh.Cells()[3].SibHfacets()[2].VertexId());
    EXPECT_EQ(-1, mesh.Cells()[3].SibHfacets()[2].CellId());
    EXPECT_EQ(-1, mesh.Cells()[3].SibHfacets()[2].FacetId());

    // Check vertex to nodes mapping.
    // Vertex Id: HalfFacet
    //    0:        <3,2>
    EXPECT_EQ(0, mesh.NodesToHalfFacets().HalfFacets()[0].VertexId());
    EXPECT_EQ(3, mesh.NodesToHalfFacets().HalfFacets()[0].CellId());
    EXPECT_EQ(2, mesh.NodesToHalfFacets().HalfFacets()[0].FacetId());
    //    1:        <1,2>
    EXPECT_EQ(1, mesh.NodesToHalfFacets().HalfFacets()[1].VertexId());
    EXPECT_EQ(1, mesh.NodesToHalfFacets().HalfFacets()[1].CellId());
    EXPECT_EQ(2, mesh.NodesToHalfFacets().HalfFacets()[1].FacetId());
    //    2:        <2,2>
    EXPECT_EQ(2, mesh.NodesToHalfFacets().HalfFacets()[2].VertexId());
    EXPECT_EQ(2, mesh.NodesToHalfFacets().HalfFacets()[2].CellId());
    EXPECT_EQ(2, mesh.NodesToHalfFacets().HalfFacets()[2].FacetId());
    //    3:        <3,2>
    EXPECT_EQ(3, mesh.NodesToHalfFacets().HalfFacets()[3].VertexId());
    EXPECT_EQ(3, mesh.NodesToHalfFacets().HalfFacets()[3].CellId());
    EXPECT_EQ(2, mesh.NodesToHalfFacets().HalfFacets()[3].FacetId());
    //    4:        <0,0>
    EXPECT_EQ(4, mesh.NodesToHalfFacets().HalfFacets()[4].VertexId());
    EXPECT_EQ(0, mesh.NodesToHalfFacets().HalfFacets()[4].CellId());
    EXPECT_EQ(0, mesh.NodesToHalfFacets().HalfFacets()[4].FacetId());
}

} // End of namespace IMP_test.
