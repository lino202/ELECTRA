#----------------

# Create the test executable.
ADD_EXECUTABLE(HalfFacetsArrayTest  half_facets_array_test.hpp  half_facets_array_test.cpp)

# Link with Google Test.
TARGET_LINK_LIBRARIES(HalfFacetsArrayTest  ${GTEST_BOTH_LIBRARIES}  ${CMAKE_THREAD_LIBS_INIT})

IF(IMP_BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(HalfFacetsArrayTest  ${PROJECT_NAME}_dynamic)
ELSE()
    TARGET_LINK_LIBRARIES(HalfFacetsArrayTest  ${PROJECT_NAME}_static)
ENDIF()

# Add the test in CTest.
ADD_TEST(HalfFacetsArrayTest  HalfFacetsArrayTest)


