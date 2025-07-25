cmake_minimum_required(VERSION 3.15)


project(DaVinci VERSION 46.12
        LANGUAGES CXX)


list(PREPEND CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake
)

# Dependencies
include(TupleToolSemiLeptonicDependencies)

# -- Subdirectories
lhcb_add_subdirectories(
    Phys/TupleToolSemiLeptonic
)

lhcb_finalize_configuration()

# Not sure if the code below is still needed
# I don't see anyhting similar in neither Phys or Analysis

if(EXISTS ${CMAKE_SOURCE_DIR}/Online/RootCnv)
  set(override_subdirs RootCnv)
endif()
# Enable Ganga integration
include(GangaTools OPTIONAL RESULT_VARIABLE GangaTools_FOUND)
if(GangaTools_FOUND)
  enable_ganga_integration()
else()
  message(WARNING "Ganga integration not available (cannot find GangaTools.cmake)")
endif()
