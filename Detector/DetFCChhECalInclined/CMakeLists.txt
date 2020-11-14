################################################################################
#Package : DetFCChhECalInclined
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCChhECalInclined)

file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetFCChhECalInclined SHARED ${sources})
target_link_libraries(DetFCChhECalInclined DD4hep::DDCore)

