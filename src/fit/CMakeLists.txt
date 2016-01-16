
# Find the required ROOT libraries
find_library(ROOT_CORE_LIB Core $ENV{ROOTSYS}/lib)
find_library(ROOT_GPAD_LIB Gpad $ENV{ROOSYS}/lib)
find_library(ROOT_HIST_LIB Hist $ENV{ROOTSYS}/lib)
find_library(ROOT_RIO_LIB RIO $ENV{ROOTSYS}/lib)
find_library(ROOT_MINUIT_LIB Minuit $ENV{ROOTSYS}/lib)
find_library(ROO_FIT RooFit $ENV{ROOTSYS}/lib)
find_library(ROO_FIT_CORE RooFitCore $ENV{ROOTSYS}/lib)

add_library(ChebyshevPoly ChebyshevPoly.cxx)

add_library(BumpHunter BumpHunter.cxx)
target_link_libraries(BumpHunter
    ${ROO_FIT}
    ${ROO_FIT_CORE}
)

add_executable(chebyshev_fit_example chebyshev_fit_example.cxx)
target_link_libraries(chebyshev_fit_example
    ${ROOT_CORE_LIB}
    ${ROOT_GPAD_LIB}
    ${ROOT_HIST_LIB}   
    ${ROOT_MINUIT_LIB}
    ChebyshevPoly 
    ${ROO_FIT}
    ${ROO_FIT_CORE}
)

add_executable(bump_hunt_evaluator bump_hunt_evaluator.cxx)
target_link_libraries(bump_hunt_evaluator
    ${ROOT_CORE_LIB}
    ${ROOT_RIO_LIB}
    ${ROOT_HIST_LIB}   
    ${ROO_FIT}
    ${ROO_FIT_CORE}
    RootFileReader
    BumpHunter
    Plotter
)

