################################################################################
# Package: DetFCChhTrackerTkLayout
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCChhTrackerTkLayout)
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/triplet DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCChhTrackerTkLayout)

file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCChhTrackerTkLayout SHARED ${sources})
target_link_libraries(DetFCChhTrackerTkLayout DD4hep::DDCore)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

