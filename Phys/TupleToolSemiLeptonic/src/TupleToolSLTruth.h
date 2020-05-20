#ifndef TupleToolSLTruth_H
#define TupleToolSLTruth_H 1

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IMCParticleTupleTool.h"  // Interface
#include "Kernel/IParticleDescendants.h"
#include "Kernel/IParticleTupleTool.h"  // Interface

#include <vector>

class IParticle2MCAssociator;

/**
 * @class TupleToolSLTruth TupleToolSLTruth.h
 * @brief Fill MC truth info if a link is present
 *
 * Uses an IParticle2MCAssociator to perform the association.
 *
 * IP2MCPAssociatorType: Implementation of IP2MCAssociator to be used. Default:
 * DaVinciSmartAssociator.
 *
 * ToolList: List of MCTupleTools to run. Default: [MCTupleToolKinematic]
 *
 * - head_TRUEID : true pid
 *
 * To add more entries, add the appropriate MCTupleTool
 * Configure the option ToolList to add MCTupleTools
 * The MCAssociation is run only once, then these tuple tools are called
 *
 * @sa DecayTreeTuple
 *
 * @author Jeremie Borel, Adlene Hicheur, Rob Lambert
 * @date   2019-11-11
 */

class TupleToolSLTruth : public TupleToolBase,
                         virtual public IParticleTupleTool {
 public:
  TupleToolSLTruth( const std::string& type, const std::string& name,
                    const IInterface* parent );

  ~TupleToolSLTruth() override = default;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                   const std::string&, Tuples::Tuple& ) override;

  StatusCode initialize() override;
  bool       isBeautyHadron( int );
  bool       isCharmHadron( int );
  bool       motherIsBeauty( LHCb::MCParticle*, const LHCb::MCParticle* );
  const LHCb::MCParticle* getMCParticle( const LHCb::Particle* );

 private:
  typedef std::vector<const LHCb::Particle*> ParticleVector;

  std::vector<IParticle2MCAssociator*> m_p2mcAssocs;
  std::vector<std::string>             m_p2mcAssocTypes;
  std::vector<std::string>             m_toolList;
  std::vector<IMCParticleTupleTool*>   m_mcTools;
  IParticleDescendants*                m_particleDescendants;
};

#endif  // TupleToolSLTruth_H
