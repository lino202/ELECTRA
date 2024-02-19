#--------------------------------------------------------------
set(LIB_NAME ${PROJECT_NAME})
string(TOLOWER ${PROJECT_NAME} LC_PROJECT_NAME)

# Find Eigen3.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find Boost.
find_package(Boost REQUIRED COMPONENTS system iostreams filesystem)
include_directories (${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

# The headers collection of the IMP library.
set(HEADERS 
    ${CMAKE_SOURCE_DIR}/src/IMP/IMP
    ${CMAKE_SOURCE_DIR}/src/IMP/Elements
    ${CMAKE_SOURCE_DIR}/src/IMP/Graphs
    ${CMAKE_SOURCE_DIR}/src/IMP/Image
    ${CMAKE_SOURCE_DIR}/src/IMP/ImageIO
    ${CMAKE_SOURCE_DIR}/src/IMP/Meshing
    ${CMAKE_SOURCE_DIR}/src/IMP/MeshIO
    ${CMAKE_SOURCE_DIR}/src/IMP/Tesselations
    ${CMAKE_SOURCE_DIR}/src/IMP/Topology
    ${CMAKE_SOURCE_DIR}/src/IMP/Utilities
    ${CMAKE_SOURCE_DIR}/src/IMP/Vectors
    ${CMAKE_SOURCE_DIR}/src/IMP/VoroIO
)

# Subdirectories of IMP.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/utilities)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/vectors)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/graphs)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/elements)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/image)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/image_io)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/mesh_io)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/meshing)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/tesselations)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/topology)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IMP/engine/voro_io)


set(${PROJECT_NAME}_LIBS 
    ${LC_PROJECT_NAME}Elements
    ${LC_PROJECT_NAME}Graphs
    ${LC_PROJECT_NAME}Image
    ${LC_PROJECT_NAME}ImageIO
    ${LC_PROJECT_NAME}Meshing 
    ${LC_PROJECT_NAME}MeshIO
    ${LC_PROJECT_NAME}Tesselations
    ${LC_PROJECT_NAME}Topology
    ${LC_PROJECT_NAME}Utilities    
    ${LC_PROJECT_NAME}Vectors
    ${LC_PROJECT_NAME}VoroIO
)

#-------- Build library --------
add_library(${LIB_NAME})
target_sources(${LIB_NAME} PRIVATE ${HEADERS})
target_link_libraries(${LIB_NAME} PUBLIC ${${PROJECT_NAME}_LIBS} Boost::filesystem)

target_include_directories(${LIB_NAME} PUBLIC
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    ${Boost_INCLUDE_DIR}
)
set_target_properties(${LIB_NAME}
    PROPERTIES
        LINKER_LANGUAGE CXX
        POSITION_INDEPENDENT_CODE 1
        DEBUG_POSTFIX "_d"
        PUBLIC_HEADER "${HEADERS}"    
        $<$<BUILD_SHARED_LIBS>:CXX_VISIBILITY_PRESET hidden>
        $<$<BUILD_SHARED_LIBS>:VISIBILITY_INLINES_HIDDEN 1>
        $<$<BUILD_SHARED_LIBS>:SOVERSION ${PROJECT_VERSION}>
        $<$<BUILD_SHARED_LIBS>:OUTPUT_NAME ${PROJECT_NAME}>
        $<$<BUILD_SHARED_LIBS>:MACOSX_RPATH ON>
        $<$<NOT:$<BUILD_SHARED_LIBS>>:ARCHIVE_OUTPUT_NAME ${PROJECT_NAME}>
)
#--------------------------------------

install(TARGETS
            ${LIB_NAME}
        EXPORT
            ${LIB_NAME}Targets
        RUNTIME
            DESTINATION ${INSTALL_BINDIR}
            COMPONENT ${PROJECT_NAME}_Runtime
        LIBRARY
            DESTINATION ${INSTALL_LIBDIR}
            COMPONENT ${PROJECT_NAME}_Runtime
        ARCHIVE
            DESTINATION ${INSTALL_LIBDIR}
            COMPONENT ${PROJECT_NAME}_Development
        PUBLIC_HEADER
            DESTINATION ${INSTALL_INCLUDEDIR}/${PROJECT_NAME}
            COMPONENT ${PROJECT_NAME}_Development
)

install(EXPORT ${LIB_NAME}Targets
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION ${INSTALL_CMAKEDIR}
    COMPONENT ${PROJECT_NAME}_Development
)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION 
        ${INSTALL_CMAKEDIR}
    PATH_VARS 
        INSTALL_INCLUDEDIR
        INSTALL_LIBDIR
)

write_basic_package_version_file(${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION 
        ${PROJECT_VERSION}
    COMPATIBILITY 
        SameMajorVersion
)

install(
    FILES 
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION 
        ${INSTALL_CMAKEDIR}
)
