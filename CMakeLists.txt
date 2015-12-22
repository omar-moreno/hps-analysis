
# Specify the minimum CMake version required
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Set the project name
project(HPS_ANALYSIS)

# Set the current version
set(HPS_ANALYSIS_VERSION_MAJOR 1)
set(HPS_ANALYSIS_VERSION_MINOR 0)

# Set the output directories of the binaries, static and shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${HPS_ANALYSIS_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${HPS_ANALYSIS_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${HPS_ANALYSIS_BINARY_DIR}/lib")

# Set additional compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Set the path to additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Include the HpsEvent headers
include_directories("$ENV{HPS_DST_PATH}/include/hps_event")

# Find ROOT
find_package(ROOT REQUIRED)

# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Include the LCIO headers
include_directories("$ENV{LCIO}/include")

# Find LibXml2
find_package(LibXml2 REQUIRED)

# include the xml2 headers
include_directories(${LIBXML2_INCLUDE_DIRS})

include_directories("$ENV{YAML_CPP_INCLUDE_PATH}")

set(UTILS_HDR_DIR "${HPS_ANALYSIS_SOURCE_DIR}/include/utils")
include_directories("${UTILS_HDR_DIR}")

# Include the analysis headers
set(ANALYSIS_HDR_DIR "${HPS_ANALYSIS_SOURCE_DIR}/include/analyses")
include_directories("${ANALYSIS_HDR_DIR}")

set(FIT_HDR_DIR "${HPS_ANALYSIS_SOURCE_DIR}/include/fit")
include_directories("${FIT_HDR_DIR}")

add_subdirectory("${HPS_ANALYSIS_SOURCE_DIR}/src/utils")

# Add the analysis subdirectory to the build
add_subdirectory("${HPS_ANALYSIS_SOURCE_DIR}/src/analyses")

# Add the fit subdirectory to the build
add_subdirectory("${HPS_ANALYSIS_SOURCE_DIR}/src/fit")

add_subdirectory("root")

