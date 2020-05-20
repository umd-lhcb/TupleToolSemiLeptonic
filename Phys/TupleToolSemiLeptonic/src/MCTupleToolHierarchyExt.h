#ifndef MCTupleToolHierarchyExt_H
#define MCTupleToolHierarchyExt_H 1

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IMCParticleTupleTool.h"  // Interface

#include "Kernel/IDaVinciAssociatorsWrapper.h"
#include "Kernel/Particle2MCLinker.h"

class MCTupleToolHierarchyExt : public TupleToolBase,
                                virtual public IMCParticleTupleTool {
 public:
  MCTupleToolHierarchyExt( const std::string& type, const std::string& name,
                           const IInterface* parent );

  virtual ~MCTupleToolHierarchyExt(){};

  StatusCode fill( const LHCb::MCParticle*, const LHCb::MCParticle*,
                   const std::string&, Tuples::Tuple& ) override;
};

#endif  // MCTupleToolHierarchyExt_H
