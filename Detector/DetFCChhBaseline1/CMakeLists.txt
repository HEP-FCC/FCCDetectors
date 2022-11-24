################################################################################
# Package: DetFCChhBaseline1
################################################################################

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Detector/DetFCChhBaseline1)




add_test(NAME ConstructDetFCChhBaseline1
         COMMAND root -q -b ".github/scripts/load_detector.C(\"Detector/DetFCChhBaseline1/compact/FCChh_DectMaster.xml\", \"fcchh_baseline.root\")"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         )
k4_set_test_env(ConstructDetFCChhBaseline1)
