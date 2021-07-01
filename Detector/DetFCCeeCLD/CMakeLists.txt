###############################################################################
# Package: DetFCCeeCLD
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCCeeCLD)


file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCCeeCLD SHARED ${sources})
target_link_libraries(DetFCCeeCLD DD4hep::DDCore)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_test(NAME ConstructDetFCCeeCLD
         COMMAND root -q -b ".github/scripts/load_detector.C(\"Detector/DetFCCeeCLD/compact/FCCee_o2_v02/FCCee_o2_v02.xml\")"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         )
k4_set_test_env(ConstructDetFCCeeCLD)
