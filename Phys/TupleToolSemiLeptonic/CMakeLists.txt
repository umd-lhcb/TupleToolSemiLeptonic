# vim: set tw=80 cc=81:
#===============================================================================
# Created    : 2018-08-21
# Maintainer : Yipeng Sun
# Package    : TupleToolSemiLeptonic
#===============================================================================

gaudi_subdir(TupleToolSemiLeptonic)
gaudi_depends_on_subdirs(Phys/DecayTreeTuple)

find_package(HepMC)
find_package(Boost)
find_package(ROOT COMPONENTS TMVA Physics)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})

gaudi_add_module(TupleToolSemiLeptonic
                 src/*.cpp
                 INCLUDE_DIRS HepMC ROOT
                 LINK_LIBRARIES HepMC DecayTreeFitter DecayTreeTupleBaseLib LoKiPhysLib ROOT)
