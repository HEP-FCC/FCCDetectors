################################################################################
# Package: DetFCCeeHCalTile
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCCeeHCalTile)

file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCCeeHCalTile SHARED ${sources})
target_link_libraries(DetFCCeeHCalTile DD4hep::DDCore)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
