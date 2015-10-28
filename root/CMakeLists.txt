
find_library(ROO_FIT RooFit $ENV{ROOTSYS}/lib)
find_library(ROO_FIT_CORE RooFitCore $ENV{ROOTSYS}/lib)

find_library(HPS_EVENT_LIBRARY HpsEvent $ENV{HPS_DST_PATH}/build/lib)

add_executable(cluster_charge_analysis cluster_charge_analysis.cxx)
target_link_libraries(cluster_charge_analysis
    RootFileReader
    ${ROO_FIT}
    ${ROO_FIT_CORE}
)

add_executable(strategy_analysis strategy_analysis.cxx)
target_link_libraries(strategy_analysis
    ${ROOT_CORE_LIBRARY}
    ${ROOT_RIO_LIBRARY}
    ${ROOT_TREE_LIBRARY}
    ${ROOT_PHYSICS_LIB}
    ${ROOT_GPAD_LIB}
    ${ROOT_HIST_LIB}
)

