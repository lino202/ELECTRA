#-----------------------------------------

# Check if Doxygen is installed in the OS.
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
    # Set input and output files.
    SET(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/imp_doxyfile.in)
    SET(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/imp_doxyfile)

    # Request to configure the Doxygen documentation file.
    CONFIGURE_FILE(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    
    # Add custom target "doc" for documentation compilation.
    ADD_CUSTOM_TARGET(doc COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          COMMENT "Generating the IMP Doxygen documentation..."
                          VERBATIM )
ELSE ()
  MESSAGE("Doxygen and Graphviz' dot need to be installed to generate the IMP Doxygen documentation.")
ENDIF ()
