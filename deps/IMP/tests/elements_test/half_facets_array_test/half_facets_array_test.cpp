#include "half_facets_array_test.hpp"

namespace IMP_test {


// Check default constructor.
TEST_F(HalfFacetsArrayTest, Ctor) {

    // Create a new HalfFacetsArray.
    IMP::HalfFacetsArray hf_array;

    // Check if is empty.
    EXPECT_EQ(hf_array.HalfFacets().size(), 0);
    EXPECT_FALSE(hf_array.IsSorted());

}


// Check Clear function.
TEST_F(HalfFacetsArrayTest, Clear) {

    // Create a new HalfFacetsArray.
    IMP::HalfFacetsArray hf_array;

    // Append some half facets.
    hf_array.Append(1,1,1);
    hf_array.Append(2,2,2);

    // Check that size is not zero.
    EXPECT_EQ(hf_array.HalfFacets().size(), 2);

    // Clear the HalfFacetsArray.
    hf_array.Clear();

    // Check if is empty.
    EXPECT_EQ(hf_array.HalfFacets().size(), 0);

}


// Check Size function.
TEST_F(HalfFacetsArrayTest, Size) {

    // Create a new HalfFacetsArray.
    IMP::HalfFacetsArray hf_array;

    // Check that size is zero.
    EXPECT_EQ(hf_array.Size(), 0);

    // Append some half facets.
    hf_array.Append(1,1,1);
    hf_array.Append(2,2,2);

    // Check that size is two.
    EXPECT_EQ(hf_array.Size(), 2);

}


// Check Append function by giving vertex, cell, facet ids.
TEST_F(HalfFacetsArrayTest, AppendWithIds) {

    // Create empty half-facets array.
    IMP::HalfFacetsArray hf_array;

    // Append first half-facet by giving its information to the array.
    hf_array.Append(1,152,4);

    // Check new size.
    EXPECT_EQ(hf_array.Size(), 1);

    // Check sorted status. Must be false.
    EXPECT_FALSE(hf_array.IsSorted());

    // Sort hf_array.
    hf_array.Sort();

    // Check sorted status. Must be true.
    EXPECT_TRUE(hf_array.IsSorted());

    // Add new half-facet in the array.
    hf_array.Append(4,125,2);

    // Check new size.
    EXPECT_EQ(hf_array.Size(), 2);

    // Check sorted status. Must be false again.
    EXPECT_FALSE(hf_array.IsSorted());
}


// Check Append function by giving half facet.
TEST_F(HalfFacetsArrayTest, AppendWithHalfFacet) {

    // Create empty half-facets array.
    IMP::HalfFacetsArray hf_array;

    // Create two half-facets.
    IMP::HalfFacet hf1(2,547,2), hf2(1520,541,3);

    // Append first half-facet by giving its information to the array.
    hf_array.Append(hf1);

    // Check new size.
    EXPECT_EQ(hf_array.Size(), 1);

    // Check sorted status. Must be false.
    EXPECT_FALSE(hf_array.IsSorted());

    // Sort hf_array.
    hf_array.Sort();

    // Check sorted status. Must be true.
    EXPECT_TRUE(hf_array.IsSorted());

    // Add new half-facet in the array.
    hf_array.Append(hf2);

    // Check new size.
    EXPECT_EQ(hf_array.Size(), 2);

    // Check sorted status. Must be false again.
    EXPECT_FALSE(hf_array.IsSorted());

}


// Check Sort function.
TEST_F(HalfFacetsArrayTest, Sort) {

    // Create a Half-facets array.
    IMP::HalfFacetsArray hf_array;

    // Append half-facets.
    hf_array.Append(2,45,23);
    hf_array.Append(5,253,1);
    hf_array.Append(35,12,5);
    hf_array.Append(2,5,25);
    hf_array.Append(2,5,18);

    // Check sorted status. Must be false.
    EXPECT_FALSE(hf_array.IsSorted());

    // Sort now the half-facets array.
    hf_array.Sort();

    // Check sorted status. Must be true.
    EXPECT_TRUE(hf_array.IsSorted());

    // Check if it is correctly sorted.
    //  2   5 18
    //  2   5 25
    //  2  45 23
    //  5 253  1
    // 35  12  5
    EXPECT_EQ(hf_array.HalfFacets()[0].VertexId(), 2);
    EXPECT_EQ(hf_array.HalfFacets()[1].CellId(), 5);
    EXPECT_EQ(hf_array.HalfFacets()[2].FacetId(), 23);
    EXPECT_EQ(hf_array.HalfFacets()[3].VertexId(), 5);
    EXPECT_EQ(hf_array.HalfFacets()[4].VertexId(), 35);

}


// Check AllAttachedToVertex function.
TEST_F(HalfFacetsArrayTest, AllAttachedToVertex) {

    // Create a Half-facets array.
    IMP::HalfFacetsArray hf_array;

    // Append half-facets.
    hf_array.Append(2,45,23);
    hf_array.Append(5,253,1);
    hf_array.Append(35,12,5);
    hf_array.Append(2,5,25);
    hf_array.Append(2,5,18);

    // Find range in array of half-facets to vertex_id = 2.
    auto range2 = hf_array.AllAttachedToVertex(2);

    // Check range.
    int diff2_1 = range2.first - hf_array.HalfFacets().begin();
    int diff2_2 = range2.second - hf_array.HalfFacets().begin();
    EXPECT_EQ(diff2_1, 0);
    EXPECT_EQ(diff2_2,3);

    // Find range in array of half-facets to vertex_id = 5.
    auto range5 = hf_array.AllAttachedToVertex(5);

    // Check range.
    int diff5_1 = range5.first - hf_array.HalfFacets().begin();
    int diff5_2 = range5.second - hf_array.HalfFacets().begin();
    EXPECT_EQ(diff5_1, 3);
    EXPECT_EQ(diff5_2,4);

    // Find range in array of half-facets to vertex_id = 10. Which doesn't exist in the array.
    auto range10 = hf_array.AllAttachedToVertex(10);

    // Check range. Since 10 doesn't belong in the array, both bound should be at the array's last element.
    int diff10_1 = range10.first - hf_array.HalfFacets().begin();
    int diff10_2 = range10.second - hf_array.HalfFacets().begin();
    int end_pos = hf_array.HalfFacets().end() - hf_array.HalfFacets().begin() - 1;
    EXPECT_EQ(diff10_1, end_pos);
    EXPECT_EQ(diff10_2, end_pos);

}


// Check FirstAttachedToVertex function.
TEST_F(HalfFacetsArrayTest, FirstAttachedToVertex) {

    // Create a Half-facets array.
    IMP::HalfFacetsArray hf_array;

    // Append half-facets.
    hf_array.Append(2,45,23);
    hf_array.Append(5,253,1);
    hf_array.Append(35,12,5);
    hf_array.Append(2,5,25);
    hf_array.Append(2,5,18);

    // Find first attached half-facet iterator to vertex_id = 2.
    auto first_attached = hf_array.FirstAttachedToVertex(2);
    int pos2 = first_attached - hf_array.HalfFacets().begin();

    // Check pos2.
    EXPECT_EQ(pos2, 0);

    // Find position in array of first attached half-facet to vertex_id = 5.
    first_attached = hf_array.FirstAttachedToVertex(5);
    int pos5 = first_attached - hf_array.HalfFacets().begin();

    // Check pos5.
    EXPECT_EQ(pos5, 3);

    // Find position in array of first attached half-facet to vertex_id = 10. Which doesn't exist in the array.
    first_attached = hf_array.FirstAttachedToVertex(10);
    int pos10 = first_attached - hf_array.HalfFacets().begin();
    int end_pos = hf_array.HalfFacets().end() - hf_array.HalfFacets().begin() -1;

    // Check pos10. Since 10 doesn't belong in the array, first attached should be the last element of the container.
    EXPECT_EQ(pos10, end_pos);
}


// Check MappedVerticesIds function.
TEST_F(HalfFacetsArrayTest, MappedVerticesIds) {

    // Create a Half-facets array.
    IMP::HalfFacetsArray hf_array;

    // Append half-facets.
    hf_array.Append(2,45,23);
    hf_array.Append(5,253,1);
    hf_array.Append(35,12,5);
    hf_array.Append(2,5,25);
    hf_array.Append(2,5,18);

    // Get the indices of the mapped vertices in the array.
    std::vector<int> mapped_verts_ids = hf_array.MappedVerticesIds();

    // Check size of the mapped_verts_ids container.
    EXPECT_EQ(mapped_verts_ids.size(), 3);

    // Check indices of mapped vertices.
    EXPECT_EQ(mapped_verts_ids[0], 2);
    EXPECT_EQ(mapped_verts_ids[1], 5);
    EXPECT_EQ(mapped_verts_ids[2], 35);

}


} // End of namespace IMP_test.
