
# Specify the minimum CMake version required
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Set the project name
project(HPS_ANALYSIS)

# Set the current version
set(HPS_DST_VERSION_MAJOR 1)
set(HPS_DST_VERSION_MINOR 0)

# Set the output directories of the binaries, static and shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${HPS_ANALYSIS_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${HPS_ANALYSIS_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${HPS_ANALYSIS_BINARY_DIR}/lib")

# Set the path to additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Find ROOT
find_package(ROOT REQUIRED)

# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Find LibXml2
find_package(LibXml2 REQUIRED)

# include the xml2 headers
include_directories(${LIBXML2_INCLUDE_DIRS})

# Include the analysis headers
set(ANALYSIS_HDR_DIR "${HPS_ANALYSIS_SOURCE_DIR}/include/analysis")
include_directories("${ANALYSIS_HDR_DIR}")

# Add the analysis subdirectory to the build
add_subdirectory("${HPS_ANALYSIS_SOURCE_DIR}/src/analysis")


