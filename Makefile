build_tool=cmake
ifeq ($(wildcard build.conf),build.conf)
include build.conf
endif

ifeq ($(build_tool),cmake)
######## CMake case ########

include $(LBUTILSROOT)/data/Makefile-cmake.mk

else
######### CMT case #########

ifeq ($(wildcard build_env.sh),build_env.sh)
CMTPROJECTPATH := $(shell bash -c ". build_env.sh ; printenv CMTPROJECTPATH")
endif

all:
%:
	@env CMTPROJECTPATH="$(CMTPROJECTPATH)" $(MAKE) -f $(LBCONFIGURATIONROOT)/data/Makefile-cmt.mk $*

endif
