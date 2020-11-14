#################################################################################
##Package : DetSegmentation
#################################################################################

file(GLOB sources src/*.cpp)
add_library(DetSegmentation SHARED ${sources})
target_link_libraries(DetSegmentation DD4hep::DDCore)
target_include_directories(DetSegmentation
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/DetSegmentation/*.h)
set_target_properties(DetSegmentation PROPERTIES PUBLIC_HEADER "${headers}")



install(TARGETS DetSegmentation
  EXPORT fccDetectorsTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/DetSegmentation"
  COMPONENT dev)

#
#include(CTest)
#gaudi_add_test(TestSegmentationPhiEta
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/phiEtaSegmentation.py)
