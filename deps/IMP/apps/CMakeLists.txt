#--------------------------------------------------------------
# Build IMP apps

add_executable(IMPmain main.cpp)
target_link_libraries(IMPmain  ${PROJECT_NAME})
target_include_directories(IMPmain PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)

add_executable(mesh_test3d mesh_test3d.cpp)
target_link_libraries(mesh_test3d  ${PROJECT_NAME})
target_include_directories(mesh_test3d PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)                  

install(TARGETS 
    IMPmain 
    mesh_test3d 
    RUNTIME 
        DESTINATION ${INSTALL_BINDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
)