################################################################################
# Package: DecayTreeTuple
################################################################################
gaudi_subdir(DecayTreeTuple)
gaudi_depends_on_subdirs(Phys/DecayTreeFitter
                         Phys/DecayTreeTupleBase
                         Phys/DecayTreeTupleDalitz
                         Phys/DecayTreeTupleJets
                         Phys/DecayTreeTupleMC
                         Phys/DecayTreeTupleMuonCalib
                         Phys/DecayTreeTupleReco
                         Phys/DecayTreeTupleTrigger
                         Phys/DecayTreeTupleTracking
                         Phys/DecayTreeTupleANNPID
                         Phys/TeslaTools
                         Phys/LoKiPhys)
find_package(HepMC)
find_package(Boost)
find_package(ROOT)
find_package(ROOT COMPONENTS TMVA Physics)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
gaudi_add_module(DecayTreeTuple
                 src/*.cpp
                 INCLUDE_DIRS HepMC ROOT
                 LINK_LIBRARIES HepMC DecayTreeFitter DecayTreeTupleBaseLib LoKiPhysLib ROOT)
gaudi_install_python_modules()
gaudi_add_test(QMTest QMTEST)
