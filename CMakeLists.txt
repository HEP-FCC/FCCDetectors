CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

project(FCCDetectors)

include(GNUInstallDirs)
include(CTest)

# CMake variable setup------------------------------------------
# Default for install prefix
# ``-DCMAKE_INSTALL_PREFIX=<location>`` when invoking CMake
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/InstallArea/ CACHE PATH
    "Install path prefix, prepended onto install directories." FORCE )
endif()

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
if(NOT CMAKE_CXX_STANDARD MATCHES "14|17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

# Find CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Dependencies--------------------------------------------------
find_package(DD4hep)


# Environment  for tests----------------------------------------
get_target_property(ddcore_lib DD4hep::DDCore LOCATION)
get_filename_component(ddcore_loc ${ddcore_lib} DIRECTORY)
function(k4_set_test_env _testname)
foreach(_subdir DetCommon DetFCCeeCLD DetFCCeeCommon DetFCCeeECalInclined DetFCCeeHCalTile DetFCCeeIDEA DetFCCeeIDEA-LAr DetFCChhBaseline1 DetFCChhCalDiscs DetFCChhECalInclined DetFCChhECalSimple DetFCChhHCalTile DetFCChhTailCatcher DetFCChhTrackerTkLayout DetSegmentation DetSensitive DetTrackerSimple DRcalo DRcomponents DRsegmentation DRsensitive )
   set( _subdirs "${CMAKE_BINARY_DIR}/Detector/${_subdir}:${_subdirs}")
endforeach()
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:${CMAKE_CURRENT_BINARY_DIR}:${ddcore_loc}:${_subdirs}:$ENV{LD_LIBRARY_PATH}")
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT "ROOT_INCLUDE_PATH=${ddcore_loc}/../include:$ENV{ROOT_INCLUDE_PATH}")
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT "DD4hepINSTALL=${ddcore_loc}/../")
endfunction()

# Subdirectories -----------------------------------------------
add_subdirectory(Detector)

# Installation and export --------------------------------------
install(FILES
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  )

install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${PROJECT_NAME}Targets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/"
  )
