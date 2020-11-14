################################################################################
# Package: DetTrackerSimple
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION Detector/DetTrackerSimple)

file(GLOB sources src/*.cpp)
add_dd4hep_plugin(DetTrackerSimple SHARED ${sources})
target_link_libraries(DetTrackerSimple DD4hep::DDCore)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

