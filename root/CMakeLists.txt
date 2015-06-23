
find_library(ROO_FIT RooFit $ENV{ROOTSYS}/lib)
find_library(ROO_FIT_CORE RooFitCore $ENV{ROOTSYS}/lib)

add_executable(cluster_charge_analysis cluster_charge_analysis.cxx)
target_link_libraries(cluster_charge_analysis
    RootFileReader
    ${ROO_FIT}
    ${ROO_FIT_CORE}
)
