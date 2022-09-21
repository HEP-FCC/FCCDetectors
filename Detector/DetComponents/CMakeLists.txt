################################################################################
# Package: DetComponents
################################################################################


file(GLOB sources src/*.cpp)
gaudi_add_module(DetComponents
                 SOURCES ${sources}
                 LINK Gaudi::GaudiAlgLib
                      k4FWCore::k4FWCore
                      DetCommon
                      EDM4HEP::edm4hep
                      DD4hep::DDCore
                      DD4hep::DDG4
                      k4FWCore::k4Interface)

install(TARGETS DetComponents
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/DetComponents"
  COMPONENT dev)


# file(GLOB sources src/*.cpp)
# add_library(DetComponents SHARED ${sources})
# target_link_libraries(DetComponents Gaudi::GaudiAlgLib
#                                     k4FWCore::k4FWCore
#                                     k4FWCore::k4Interface
#                                     EDM4HEP::edm4hep
#                                     DD4hep::DDCore
#                                     DD4hep::DDG4
#                                     DetCommon)
# target_include_directories(DetComponents
#     PUBLIC
#         $<INSTALL_INTERFACE:include>
#         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
# )
# 
# file(GLOB headers include/DetComponents/*.h)
# set_target_properties(DetComponents PROPERTIES PUBLIC_HEADER "${headers}")
# dd4hep_generate_rootmap(DetComponents)
# 
# 
# 
# install(TARGETS DetComponents
#   EXPORT ${PROJECT_NAME}Targets
#   RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
#   LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
#   PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/DetComponents"
#   COMPONENT dev)
# 
# 
# set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})




#include(CTest)
#gaudi_add_test(GeantFullSimGdml
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/geant_fullsim_gdml.py)
#gaudi_add_test(GeantFullSimHCal
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/geant_fullsim_hcal.py)
#gaudi_add_test(GeantFullSimMoreEvents
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Sim/SimG4Components/tests/
#               FRAMEWORK tests/options/geant_fullsim_moreEvents.py)
#gaudi_add_test(GeantFullSimMoreEventsCheckNumEvents
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Sim/SimG4Components/tests/
#               COMMAND python ./scripts/geant_fullsim_moreEvents_checkNumEvents.py
#               DEPENDS GeantFullSimMoreEvents)
#gaudi_add_test(GeantFullSimMoreEventsCheckNumParticles
#               ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Sim/SimG4Components/tests/
#               COMMAND python ./scripts/geant_fullsim_moreEvents_checkNumParticles.py
#               DEPENDS GeantFullSimMoreEvents)
#gaudi_add_test(GeantFastSimSimpleSmearing
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Sim/SimG4Components/tests/
#               FRAMEWORK tests/options/geant_fastsim_simple.py)
#add_test(      NAME  SimG4Components.GeantFastSimSimpleSmearingCheckSigma
#               COMMAND python ./scripts/geant_fastsim_simple_checkSigma.py
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Sim/SimG4Components/tests/ )
#set_tests_properties(SimG4Components.GeantFastSimSimpleSmearingCheckSigma PROPERTIES DEPENDS SimG4Components.GeantFastSimSimpleSmearing)
#gaudi_add_test(GeantFastSimSimpleSmearingCheckNumParticles
#               ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Sim/SimG4Components/tests/
#               COMMAND python ./scripts/geant_fastsim_checkNumParticles.py
#               DEPENDS GeantFastSimSimpleSmearing)
#
