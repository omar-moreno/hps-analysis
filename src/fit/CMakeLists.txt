
# Find the required ROOT libraries
find_library(ROOT_CORE_LIB Core $ENV{ROOTSYS}/lib)
find_library(ROOT_GPAD_LIB Gpad $ENV{ROOSYS}/lib)
find_library(ROOT_HIST_LIB Hist $ENV{ROOTSYS}/lib)
find_library(ROOT_RIO_LIB RIO $ENV{ROOTSYS}/lib)
find_library(ROOT_MINUIT_LIB Minuit $ENV{ROOTSYS}/lib)
find_library(ROO_FIT RooFit $ENV{ROOTSYS}/lib)
find_library(ROO_FIT_CORE RooFitCore $ENV{ROOTSYS}/lib)
find_library(ROOT_MATH_LIB MathCore $ENV{ROOTSYS}/lib)
find_library(ROOT_TREE_LIB Tree $ENV{ROOTSYS}/lib)

add_library(ChebyshevPoly ChebyshevPoly.cxx)

add_library(HpsFitResult HpsFitResult.cxx
    ${ROO_FIT}
    ${ROO_FIT_CORE}
)

add_library(BumpHunter BumpHunter.cxx)
target_link_libraries(BumpHunter
    ${ROOT_GPAD_LIB}
    ${ROO_FIT}
    ${ROO_FIT_CORE}
    ${ROOT_MATH_LIB}
    HpsFitResult
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
    FlatTupleMaker
    HpsFitResult
)

add_executable(set_upper_limits set_upper_limits.cxx)
target_link_libraries(set_upper_limits
    ${ROOT_CORE_LIB}
    ${ROOT_RIO_LIB}
    ${ROOT_HIST_LIB} 
    ${ROOT_TREE_LIB} 
    ${ROO_FIT}
    ${ROO_FIT_CORE}
    BumpHunter
    FlatTupleMaker
    HpsFitResult
)
