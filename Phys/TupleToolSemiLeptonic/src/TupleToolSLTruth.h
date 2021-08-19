// $Id: TupleToolSLTruth.h,v 1.9 2010-01-26 15:39:26 rlambert Exp $
#ifndef TUPLETOOLMCTRUTH_H
#define TUPLETOOLMCTRUTH_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IMCParticleTupleTool.h"  // Interface
#include "Kernel/IParticleDescendants.h"
#include "Kernel/IParticleTupleTool.h"  // Interface
// struct Particle2MCLinker;
#include <vector>

class IParticle2MCAssociator;

/** @class TupleToolSLTruth TupleToolSLTruth.h
 *
 * \brief Fill MC truth info if a link is present
 *
 * Uses an IParticle2MCAssociator to perform the association.
 * <b> Properties: </b>
 *
 * IP2MCPAssociatorType: Implementation of IP2MCAssociator to be used. Default:
 DaVinciSmartAssociator.
 *
 * ToolList: List of MCTupleTools to run. Default: [MCTupleToolKinematic]
 *
 * - head_TRUEID : true pid
 *

 * To add more entries, add the appropriate MCTupleTool

 * Configure the option ToolList to add MCTupleTools

 * The MCAssociation is run only once, then these tuple tools are called


 * \sa DecayTreeTuple
 *
 *  @author Jeremie Borel
 *  @date   2007-11-07
 *  2008-09-23 Adlene Hicheur - Added true angles information for P2VV
 *  2009-06-03 Rob Lambert - Major Changes
 */

class TupleToolSLTruth : public TupleToolBase,
                         virtual public IParticleTupleTool {
 public:
  /// Standard constructor
  TupleToolSLTruth( const std::string& type, const std::string& name,
                    const IInterface* parent );

  virtual ~TupleToolSLTruth(){};  ///< Destructor

  StatusCode getMCTruePi0s( const LHCb::MCParticle* MCP,
                            std::vector<double>*    MCTruePi0PX_,
                            std::vector<double>*    MCTruePi0PY_,
                            std::vector<double>*    MCTruePi0PZ_,
                            std::vector<double>*    MCTruePi0E_,
                            std::vector<int>*       MCTruePi0MID_ );

  StatusCode getMCTrueGammas( const LHCb::MCParticle* MCP,
                              std::vector<double>*    MCTrueGammaPX_,
                              std::vector<double>*    MCTrueGammaPY_,
                              std::vector<double>*    MCTrueGammaPZ_,
                              std::vector<double>*    MCTrueGammaE_,
                              std::vector<int>*       MCTrueGammaMID_,
                              std::vector<int>*       MCTrueGammaIsRad_ );

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                   const std::string&, Tuples::Tuple& ) override;

  StatusCode initialize() override;
  bool       isBeautyHadron( int );
  bool       isCharmHadron( int );
  bool       motherIsBeauty( LHCb::MCParticle*, const LHCb::MCParticle* );
  const LHCb::MCParticle* getMCParticle( const LHCb::Particle* );

 private:
  typedef std::vector<const LHCb::Particle*> ParticleVector;
  std::vector<IParticle2MCAssociator*>       m_p2mcAssocs;
  std::vector<std::string>                   m_p2mcAssocTypes;
  std::vector<std::string>
      m_toolList;  ///< names of all MCTupleTools, set by the option ToolList

  std::vector<IMCParticleTupleTool*> m_mcTools;  ///< vector of MCTools to fill
  IParticleDescendants*              m_particleDescendants;
};

#endif  // TUPLETOOLMCTRUTH_H
