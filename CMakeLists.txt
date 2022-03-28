CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

project(FCCDetectors)




include(GNUInstallDirs)
include(CTest)

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

message(${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
message(${CMAKE_MODULE_PATH})


#---------------------------------------------------------------
find_package(DD4hep)
#---------------------------------------------------------------

get_target_property(ddcore_lib DD4hep::DDCore LOCATION)
get_filename_component(ddcore_loc ${ddcore_lib} DIRECTORY)
function(k4_set_test_env _testname)
foreach(_subdir DetCommon DetFCCeeCLD DetFCCeeCommon DetFCCeeECalInclined DetFCCeeHCalTile DetFCCeeIDEA DetFCCeeIDEA-LAr DetFCChhBaseline1 DetFCChhCalDiscs DetFCChhECalInclined DetFCChhECalSimple DetFCChhHCalTile DetFCChhTailCatcher DetFCChhTrackerTkLayout DetSegmentation DetSensitive DetTrackerSimple)
   set( _subdirs "${CMAKE_BINARY_DIR}/Detector/${_subdir}:${_subdirs}")
endforeach()
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:${CMAKE_CURRENT_BINARY_DIR}:${ddcore_loc}:${_subdirs}:$ENV{LD_LIBRARY_PATH}")
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT "ROOT_INCLUDE_PATH=${ddcore_loc}/../include:$ENV{ROOT_INCLUDE_PATH}")
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT "DD4hepINSTALL=${ddcore_loc}/../")

endfunction()



add_subdirectory(Detector)

# - install and export
install(FILES
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  )






install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${PROJECT_NAME}Targets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/"
  )

#gaudi_install(CMAKE cmake/${PROJECT_NAME}Config.cmake)

