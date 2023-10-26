################################################################################
#Package : DetFCCeeECalInclined
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCCeeECalInclined)

file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCCeeECalInclined SHARED ${sources})
target_link_libraries(DetFCCeeECalInclined DD4hep::DDCore)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
