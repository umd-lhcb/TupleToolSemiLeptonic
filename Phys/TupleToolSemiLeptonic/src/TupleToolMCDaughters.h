// $Id: TupleToolMCDaughters.h,v 1.9 2010-01-26 15:39:26 rlambert Exp $
#ifndef TupleToolMCDaughters_H
#define TupleToolMCDaughters_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IMCParticleTupleTool.h"  // Interface
#include "Kernel/IParticleTupleTool.h"    // Interface

//struct Particle2MCLinker;
#include <vector>

class IParticle2MCAssociator;

class TupleToolMCDaughters : public TupleToolBase, virtual public IParticleTupleTool
{

public:

  /// Standard constructor
  TupleToolMCDaughters( const std::string& type,
                    const std::string& name,
                    const IInterface* parent);

  virtual ~TupleToolMCDaughters(){}; ///< Destructor

  virtual StatusCode fill( const LHCb::Particle*
                           , const LHCb::Particle*
                           , const std::string&
                           , Tuples::Tuple& ) override;

  virtual StatusCode initialize() override;

private:

  std::vector<IParticle2MCAssociator*> m_p2mcAssocs;
  std::vector<std::string> m_p2mcAssocTypes;
  bool m_Mother;
  unsigned m_nMin;
  std::vector<std::string> m_toolList; ///< names of all MCTupleTools, set by the option ToolList

  std::vector<IMCParticleTupleTool*> m_mcTools; ///<vector of MCTools to fill

};

#endif // TupleToolMCDaughters_H
