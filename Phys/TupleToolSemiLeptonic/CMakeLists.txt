# vim: set tw=80 cc=81:
#===============================================================================
# Created    : 2019-08-21
# Maintainer : Yipeng Sun
# Package    : TupleToolSemiLeptonic
#===============================================================================

gaudi_add_module(TupleToolSemiLeptonic
    SOURCES
        src/MCTupleToolHierarchyExt.cpp
        src/TupleToolApplyIsolation.cpp
        src/TupleToolApplyIsolationVetoDst.cpp
        src/TupleToolDocas.cpp
        src/TupleToolMCDaughters.cpp
        src/TupleToolSLTools.cpp
        src/TupleToolSLTruth.cpp
        src/TupleToolTagDiscardDstMu.cpp
        src/TupleToolTauMuDiscrVars.cpp
    LINK
        Analysis::DecayTreeTupleBaseLib
        Phys::DecayTreeFitter
        Phys::LoKiPhysLib
        ROOT::Core
        ROOT::Physics
        ROOT::RIO
        ROOT::TMVA
)
