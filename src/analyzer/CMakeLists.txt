
find_library(HPS_EVENT_LIBRARY HpsEvent $ENV{HPS_DST_PATH}/build/lib)

add_executable(hps_analysis hps_analysis.cxx)
target_link_libraries(hps_analysis 
    ${ROOT_CORE_LIBRARY}
    ${ROOT_RIO_LIBRARY}
    ${ROOT_TREE_LIBRARY}
    ${LIBXML2_LIBRARIES}
    ${HPS_EVENT_LIBRARY}
    TrackAnalysis
    PairsAnalysis
    TagProbeAnalysis
    TrackClusterMatchingEfficiencyAnalysis
    SharedHitAnalysis
    TimingAnalysis
    MollerAnalysis
    MuonAnalysis
    GblTrackAnalysis
    TridentAnalysis
)

# Find the necessary LCIO libraries
find_library(LCIO_LIB lcio $ENV{LCIO}/lib)

add_executable(lcio_analyzer lcio_analyzer.cxx)
target_link_libraries(lcio_analyzer ${LCIO_LIB})

