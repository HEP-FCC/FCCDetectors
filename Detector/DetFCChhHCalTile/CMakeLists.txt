################################################################################
# Package: DetFCChhHCalTile
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCChhHCalTile)


file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCChhHCalTile SHARED ${sources})
target_link_libraries(DetFCChhHCalTile DD4hep::DDCore)


set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


add_test(NAME ConstructDetFCChhHCalTile
         COMMAND root -q -b ".github/scripts/load_detector.C(\"Detector/DetFCChhHCalTile/compact/standalone.xml\")"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         )
k4_set_test_env(ConstructDetFCChhHCalTile)
