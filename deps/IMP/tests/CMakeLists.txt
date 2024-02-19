#----------------------------------------------
NAME test_shared
COMMAND $<TARGET_FILE:hello-world_wDSO>
)

add_test(
NAME test_static
COMMAND $<TARGET_FILE:hello-world_wAR>
)

add_test(
NAME use-target_configure
COMMAND
  ${CMAKE_COMMAND} -H${CMAKE_CURRENT_LIST_DIR}/use_target
                   -B${CMAKE_CURRENT_BINARY_DIR}/build_use-target
                   -G${CMAKE_GENERATOR}
                   -Dmessage_DIR=${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKEDIR}
                   -DCMAKE_BUILD_TYPE=$<CONFIGURATION>
)
set_tests_properties(use-target_configure
PROPERTIES
  FIXTURES_SETUP use-target
)

add_test(
NAME use-target_build
COMMAND
  ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/build_use-target
                   --config $<CONFIGURATION>
)
set_tests_properties(use-target_build
PROPERTIES
  FIXTURES_REQUIRED use-target
)

set(_test_target)
if(MSVC)
set(_test_target "RUN_TESTS")
else()
set(_test_target "test")
endif()
add_test(
NAME use-target_test
COMMAND
  ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/build_use-target
                   --target ${_test_target}
                   --config $<CONFIGURATION>
)
set_tests_properties(use-target_test
PROPERTIES
  FIXTURES_REQUIRED use-target
)
unset(_test_target)

add_test(
NAME use-target_cleanup
COMMAND
  ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/build_use-target
)
set_tests_properties(use-target_cleanup
PROPERTIES
  FIXTURES_CLEANUP use-target
)



# # Add definitions for compilation in Apple Architecture.
# IF (APPLE)
#     ADD_DEFINITIONS(-DGTEST_USE_OWN_TR1_TUPLE)
#     ADD_DEFINITIONS(-D__GLIBCXX__)
# ENDIF (APPLE)

# # Find Google Test dependency.
# FIND_PACKAGE(GTest  REQUIRED)
# INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})

# # Copy image testing data to binary directory.
# FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_double.mhd
#           ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_double.raw
#           ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_float.mhd
#           ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_float.raw
#           ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_int.mhd
#           ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_int.raw
#           ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_uint.mhd
#           ${CMAKE_CURRENT_SOURCE_DIR}/image_test/test_data/test_img_uint.raw
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/image_test/test_data)

# # Add testing directories.
# ADD_SUBDIRECTORY(vectors_test/vec_test)
# ADD_SUBDIRECTORY(image_test/voxel_image_test)
# ADD_SUBDIRECTORY(elements_test/cell_test)
# ADD_SUBDIRECTORY(elements_test/half_facet_test)
# ADD_SUBDIRECTORY(elements_test/half_facets_array_test)
# ADD_SUBDIRECTORY(meshes_test/mesh_test)
