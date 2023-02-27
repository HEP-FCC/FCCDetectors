###############################################################################
# Package: DetFCCeeIDEA-LAr
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCCeeIDEA-LAr)

add_test(NAME ConstructDetFFCCeeIDEALAr
         COMMAND root -q -b ".github/scripts/load_detector.C(\"Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster.xml\", \"${CMAKE_BINARY_DIR}/fccee_nlq.root\")"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         )
k4_set_test_env(ConstructDetFFCCeeIDEALAr)


