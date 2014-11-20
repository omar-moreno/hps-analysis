
# Find the necessary LCIO libraries
find_library(LCIO_LIB lcio $ENV{LCIO}/lib)

add_executable(lcio_analyzer lcio_analyzer.cxx)
target_link_libraries(lcio_analyzer ${LCIO_LIB})

