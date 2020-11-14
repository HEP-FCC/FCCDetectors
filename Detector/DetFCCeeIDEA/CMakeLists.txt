###############################################################################
# Package: DetFCCeeIDEA
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCCeeIDEA)


file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCCeeIDEA SHARED ${sources})
target_link_libraries(DetFCCeeIDEA DD4hep::DDCore DetSegmentation)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
