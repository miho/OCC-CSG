CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

PROJECT(occ-csg-prj)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

INCLUDE(cmake/CreateExecutable.cmake)
INCLUDE(cmake/CreateLibrary.cmake)

# Find opencascade libraries

find_package(OpenCASCADE COMPONENTS FoundationClasses REQUIRED)
find_package(OpenCASCADE COMPONENTS ModelingAlgorithms REQUIRED)
find_package(OpenCASCADE COMPONENTS ModelingData REQUIRED)
find_package(OpenCASCADE COMPONENTS DataExchange REQUIRED)
find_package(OpenCASCADE COMPONENTS Visualization REQUIRED)
include_directories(${OpenCASCADE_INCLUDE_DIR})
MESSAGE("OCC-Include-Dir: ${OpenCASCADE_INCLUDE_DIR}")
link_directories(${OpenCASCADE_LIBRARY_DIR})
MESSAGE("OCC-Lib-Dir-Dir: ${OpenCASCADE_LIBRARY_DIR}")


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/dependencies/freetype/include")

# Find openmp libraries

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# Find freetype libraries on Linux (don't know why the one provided by the static occt dep doesn't work on Linux)
IF(UNIX AND NOT APPLE)
    find_library(FreeType REQUIRED)
#    include_directories(${FREETYPE_INCLUDE_DIRS})
#    link_directories(${FREETYPE_LIBRARY_DIR})
ENDIF()

set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

ADD_SUBDIRECTORY(src)
