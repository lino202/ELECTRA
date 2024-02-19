#----------------

# Set the test's name.
SET(TestName  MeshTest)

# Create the test executable.
ADD_EXECUTABLE(${TestName}  mesh_test.hpp  mesh_test.cpp)

# Link with Google Test.
TARGET_LINK_LIBRARIES(${TestName}  ${GTEST_BOTH_LIBRARIES}  ${CMAKE_THREAD_LIBS_INIT})

IF(IMP_BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(${TestName}  ${PROJECT_NAME}_dynamic)
ELSE()
    TARGET_LINK_LIBRARIES(${TestName}  ${PROJECT_NAME}_static)
ENDIF()

# Add the test in CTest.
ADD_TEST(${TestName}  ${TestName})


