
find_library(HPS_EVENT_LIBRARY HpsEvent $ENV{HPS_DST_PATH}/build/lib)

add_library(TrackAnalysis TrackAnalysis.cxx)
target_link_libraries(TrackAnalysis
    ${ROOT_CORE_LIBRARY}
    ${ROOT_GPAD_LIBRARY}
    ${ROOT_HIST_LIBRARY}
    ${HPS_EVENT_LIBRARY}
)