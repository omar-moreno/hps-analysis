# Find the required ROOT libraries
find_library(ROOT_TREE_LIB Tree $ENV{ROOTSYS}/lib)
find_library(ROOT_CORE_LIB Core $ENV{ROOTSYS}/lib)
find_library(ROOT_RIO_LIB RIO $ENV{ROOTSYS}/lib)
find_library(ROOT_PHYSICS_LIB Physics $ENV{ROOTSYS}/lib)
find_library(ROOT_HIST_LIB Hist $ENV{ROOTSYS}/lib)
find_library(ROOT_GPAD_LIB Gpad $ENV{ROOTSYS}/lib)
find_library(ROOT_GRAF_LIB Graf $ENV{ROOTSYS}/lib)

find_library(HPS_EVENT_LIBRARY HpsEvent $ENV{HPS_DST_PATH}/build/lib)

find_library(YAML_CPP_LIB yaml-cpp $ENV{YAML_CPP_LIB_PATH})

add_library(TrackType TrackType.cxx)
target_link_libraries(TrackType
    ${HPS_EVENT_LIBRARY}    
)

add_library(TrackUtils TrackUtils.cxx)
target_link_libraries(TrackUtils
    ${HPS_EVENT_LIBRARY}
    TrackType
)

add_library(AnalysisUtils AnalysisUtils.cxx)
target_link_libraries(AnalysisUtils
    ${HPS_EVENT_LIBRARY}
    TrackExtrapolator
)

add_library(Plotter Plotter.cxx) 
target_link_libraries(Plotter
    ${ROOT_CORE_LIB}
    ${ROOT_HIST_LIB}
)

add_library(PlottingUtils PlottingUtils.cxx)
target_link_libraries(PlottingUtils 
    ${ROOT_CORE_LIB}
)

add_library(TrackExtrapolator TrackExtrapolator.cxx)
target_link_libraries(TrackExtrapolator
    ${HPS_EVENT_LIBRARY}
)

add_library(EcalUtils EcalUtils.cxx)
target_link_libraries(EcalUtils
    ${HPS_EVENT_LIBRARY}
    Plotter
)


add_library(TrackClusterMatcher TrackClusterMatcher.cxx)
target_link_libraries(TrackClusterMatcher
    ${HPS_EVENT_LIBRARY}
    ${ROOT_CORE_LIB}
    ${ROOT_HIST_LIB}
    TrackExtrapolator
    Plotter
    TrackType
    TrackUtils
    AnalysisUtils
)

add_library(RootFileReader RootFileReader.cxx)
target_link_libraries(RootFileReader 
    ${ROOT_CORE_LIB}
    ${ROOT_GPAD_LIB}
    ${ROOT_TREE_LIB}
    ${ROOT_RIO_LIB}
    ${ROOT_PHYSICS_LIB}
    ${ROOT_HIST_LIB}
)

add_library(ComparePlots ComparePlots.cxx)
target_link_libraries(ComparePlots 
    ${ROOT_CORE_LIB}
    ${ROOT_GPAD_LIB}
    ${ROOT_TREE_LIB}
    ${ROOT_RIO_LIB}
    ${ROOT_PHYSICS_LIB}
    ${ROOT_HIST_LIB}
    ${ROOT_GRAF_LIB}
    RootFileReader
)

add_executable(compare_plots compare_plots.cxx)
target_link_libraries(compare_plots 
   ${YAML_CPP_LIB}
    ComparePlots    
)

add_executable(print_plots print_plots.cxx)
target_link_libraries(print_plots
    ${ROOT_CORE_LIB}
    ${ROOT_HIST_LIB}
    RootFileReader
    PlottingUtils
)

add_library(RooFitter RooFitter.cxx)
target_link_libraries(RooFitter
    ${ROO_FIT}
    ${ROO_FIT_CORE}
    ${ROOT_HIST_LIB}
)
