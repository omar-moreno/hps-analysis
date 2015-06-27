# Find the required ROOT libraries
find_library(ROOT_TREE_LIB Tree $ENV{ROOTSYS}/lib)
find_library(ROOT_CORE_LIB Core $ENV{ROOTSYS}/lib)
find_library(ROOT_RIO_LIB RIO $ENV{ROOTSYS}/lib)
find_library(ROOT_PHYSICS_LIB Physics $ENV{ROOTSYS}/lib)
find_library(ROOT_HIST_LIB Hist $ENV{ROOTSYS}/lib)
find_library(ROOT_GPAD_LIB Gpad $ENV{ROOTSYS}/lib)
find_library(ROOT_GRAF_LIB Graf $ENV{ROOTSYS}/lib)

find_library(HPS_EVENT_LIBRARY HpsEvent $ENV{HPS_DST_PATH}/build/lib)

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
    ComparePlots    
)
