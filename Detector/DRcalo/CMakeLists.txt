################################################################################
# Package: Dual Readout Calorimeter
################################################################################

file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
)

dd4hep_add_plugin(ddDRcalo SOURCES src/*.cpp USES
  DD4hep::DDCore
  DD4hep::DDCond
  ROOT::Core
  ROOT::Geom
  ROOT::GenVector
  ROOT::MathCore
  DRsegmentation
)

target_include_directories(ddDRcalo PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(ddDRcalo PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS ddDRcalo
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DRcalo)

# This test scan the material along the X (or Z) axis, 
# then extract the second-to-last line, which contains a summary of the scan, 
# then fail if the argument 12, which correspond to the integrated X0 is less than the X0 of air

add_test(NAME ScanMaterialDRcaloX
        COMMAND bash ".github/scripts/materialScanAxis.sh" "Detector/DRcalo/compact/DRcalo.xml" "x" "noair"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         )
add_test(NAME ScanMaterialDRcaloZ
         COMMAND bash ".github/scripts/materialScanAxis.sh" "Detector/DRcalo/compact/DRcalo.xml" "z" "air"
          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
          )

k4_set_test_env(ScanMaterialDRcaloX)
k4_set_test_env(ScanMaterialDRcaloZ)
