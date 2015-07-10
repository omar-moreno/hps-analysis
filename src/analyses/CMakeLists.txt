
find_library(HPS_EVENT_LIBRARY HpsEvent $ENV{HPS_DST_PATH}/build/lib)

add_library(TrackAnalysis TrackAnalysis.cxx)
target_link_libraries(TrackAnalysis
    ${ROOT_CORE_LIBRARY}
    ${ROOT_RIO_LIBRARY}
    ${ROOT_GPAD_LIBRARY}
    ${ROOT_HIST_LIBRARY}
    ${HPS_EVENT_LIBRARY}
    Plotter
    TrackExtrapolator
)

add_library(PairsAnalysis PairsAnalysis.cxx)
target_link_libraries(PairsAnalysis
    ${ROOT_CORE_LIBRARY}
    ${ROOT_RIO_LIBRARY}
    ${ROOT_GPAD_LIBRARY}
    ${ROOT_HIST_LIBRARY}
    ${HPS_EVENT_LIBRARY}
    Plotter
)

add_library(TagProbeAnalysis TagProbeAnalysis.cxx)
target_link_libraries(TagProbeAnalysis
    ${HPS_EVENT_LIBRARY}
    TrackExtrapolator
    Plotter
)

add_library(SimpleTrackingEfficiencyAnalysis SimpleTrackingEfficiencyAnalysis.cxx)
target_link_libraries(SimpleTrackingEfficiencyAnalysis
    ${HPS_EVENT_LIBRARY}
    TrackExtrapolator
    Plotter
)

add_library(SharedHitAnalysis SharedHitAnalysis.cxx)
target_link_libraries(SharedHitAnalysis
    ${HPS_EVENT_LIBRARY}
    Plotter
)

add_library(TimingAnalysis TimingAnalysis.cxx)
target_link_libraries(SimpleTrackingEfficiencyAnalysis
    ${HPS_EVENT_LIBRARY}
    ${ROOT_CORE_LIBRARY}
    ${ROOT_HIST_LIBRARY}
    ${ROOT_GPAD_LIBRARY}
    Plotter
)
