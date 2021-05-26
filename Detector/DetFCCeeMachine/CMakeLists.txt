###############################################################################
# Package: DetFCCeeMachine
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCCeeMachine)


file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCCeeMachine SHARED ${sources})
target_link_libraries(DetFCCeeMachine DD4hep::DDCore DetSegmentation)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


add_test(NAME ConstructDetFFCCeeMachine
         COMMAND root -q -b ".github/scripts/load_detector.C(\"Detector/DetFCCeeMachine/compact/FCCee_DectMaster.xml\")"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         )
k4_set_test_env(ConstructDetFFCCeeMachine)


