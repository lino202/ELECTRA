# - Find Voro++
#
# Find the native Voro++ includes and add_library
#
# Voro++_FOUND - True if Voro++ found.
# Voro++_INCLUDE_DIR - Where to find voro++.hh, etc.
# Voro++_INCLUDE_DIRS - Set if found voro++.hh, etc.
# Voro++_LIBRARIES - List of libraries when using Voro++.
#

include("FindPackageHandleStandardArgs")

if (MSVC)
    set(Voro++_DIRS "C:/Program Files (x86)" CACHE PATH "Directory containing Voro++ directories")
else ()
    set(Voro++_DIRS "/usr/local" CACHE PATH "Directory containing Voro++ directories")
endif ()

# Find header files.
find_path(Voro++_INCLUDE_DIR "voro++.hh" PATHS ${Voro++_DIRS}/include PATH_SUFFIXES "voro++")

#Find libraries.
find_library(Voro++_LIBRARY NAMES "voro++" PATHS ${Voro++_DIRS}/lib PATH_SUFFIXES "voro++")

# handle the QUIETLY and REQUIRED arguments and set Voro++_FOUND to TRUE
# if all listed variables are TRUE.

FIND_PACKAGE_HANDLE_STANDARD_ARGS("Voro++" DEFAULT_MSG Voro++_INCLUDE_DIR Voro++_LIBRARY)

if(Voro++_FOUND)
    set(Voro++_INCLUDE_DIRS  ${Voro++_INCLUDE_DIR})
    set(Voro++_LIBRARIES  ${Voro++_LIBRARY})
else()
    message(FATAL_ERROR "Could not locate Voro++...")
endif()
