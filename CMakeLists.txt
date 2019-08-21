CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

#---------------------------------------------------------------
# Load macros and functions for Gaudi-based projects
find_package(GaudiProject)
#---------------------------------------------------------------

if(EXISTS ${CMAKE_SOURCE_DIR}/Online/RootCnv)
  set(override_subdirs RootCnv)
endif()

# Declare project name and version
# Syntax:
#    gaudi_project(this_project this_version
#                  USE dep_project_1 version_1 [project_2 version_2 ...]
#                  [DATA pkg1 [VERSION vers1] [pkg2 ...])
gaudi_project(DaVinciDev  v45r1
              USE DaVinci v45r1)

# Enable Ganga integration
include(GangaTools OPTIONAL RESULT_VARIABLE GangaTools_FOUND)
if(GangaTools_FOUND)
  enable_ganga_integration()
else()
  message(WARNING "Ganga integration not available (cannot find GangaTools.cmake)")
endif()
