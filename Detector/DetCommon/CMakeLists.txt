################################################################################
# Package: DetCommon
################################################################################

file(GLOB sources src/*.cpp)
add_library(DetCommon SHARED ${sources})
target_link_libraries(DetCommon DD4hep::DDCore DD4hep::DDG4 DetSegmentation)
target_include_directories(DetCommon
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/DetCommon/*.h)
set_target_properties(DetCommon PROPERTIES PUBLIC_HEADER "${headers}")



install(TARGETS DetCommon
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/DetCommon"
  COMPONENT dev)


install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetCommon)

dd4hep_generate_rootmap(DetCommon)

#include(CTest)
#gaudi_add_test(DumpSimpleBox
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/dumpSimpleBox.py)
