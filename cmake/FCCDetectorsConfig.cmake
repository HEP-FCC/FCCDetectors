
#include(CMakeFindDependencyMacro)
#find_dependency(DD4hep REQUIRED)

# - Include the targets file to create the imported targets that a client can
# link to (libraries) or execute (programs)
include("${CMAKE_CURRENT_LIST_DIR}/FCCDetectorsTargets.cmake")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FCCDetectors  DEFAULT_MSG CMAKE_CURRENT_LIST_FILE)
