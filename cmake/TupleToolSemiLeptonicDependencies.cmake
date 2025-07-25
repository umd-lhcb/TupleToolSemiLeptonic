if(NOT COMMAND lhcb_find_package)
  # Look for LHCb find_package wrapper
  find_file(LHCbFindPackage_FILE LHCbFindPackage.cmake)
  if(LHCbFindPackage_FILE)
      include(${LHCbFindPackage_FILE})
  else()
      # if not found, use the standard find_package
      macro(lhcb_find_package)
          find_package(${ARGV})
      endmacro()
  endif()
endif()

# -- Public dependencies
lhcb_find_package(DaVinci REQUIRED)

find_package(ROOT REQUIRED
    Core
    Physics
    RIO
    TMVA
)
