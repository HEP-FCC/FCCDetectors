################################################################################
#Package : DetFCChhCalDiscs
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCChhCalDiscs)

file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCChhCalDiscs SHARED ${sources})
target_link_libraries(DetFCChhCalDiscs DD4hep::DDCore)
