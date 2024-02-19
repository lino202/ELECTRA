#----------------

# Create the test executable.
ADD_EXECUTABLE(VecTest vec_test.hpp vec_test.cpp)

# Link with Google Test.
TARGET_LINK_LIBRARIES(VecTest ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

IF(IMP_BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(VecTest  ${PROJECT_NAME}_dynamic)
ELSE()
    TARGET_LINK_LIBRARIES(VecTest  ${PROJECT_NAME}_static)
ENDIF()

# Add the test in CTest.
ADD_TEST(VecTest  VecTest)


