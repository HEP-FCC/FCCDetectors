###############################################################################
# Package: DetFCCeeCommon
################################################################################
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCCeeCommon)

file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCCeeCommon SHARED ${sources})
target_link_libraries(DetFCCeeCommon DD4hep::DDCore)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
