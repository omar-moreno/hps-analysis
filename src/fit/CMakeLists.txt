
# Find the required ROOT libraries
find_library(ROOT_CORE_LIB Core $ENV{ROOTSYS}/lib)
find_library(ROOT_GPAD_LIB Gpad $ENV{ROOSYS}/lib)
find_library(ROOT_HIST_LIB Hist $ENV{ROOTSYS}/lib)

add_executable(chebyshev_fit_example chebyshev_fit_example.cxx)
target_link_libraries(chebyshev_fit_example
    ${ROOT_CORE_LIB}
    ${ROOT_GPAD_LIB}
    ${ROOT_HIST_LIB}    
)
