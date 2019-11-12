#ifndef TupleToolApplyIsolation_H
#define TupleToolApplyIsolation_H 1

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticle2MCAssociator.h"
#include "Kernel/IParticleTupleTool.h"  // Interface
#include "Kernel/Particle2MCLinker.h"
#include "TrackInterfaces/ITrackVertexer.h"

#include "TMVA/Reader.h"
#include "TString.h"

class IDVAlgorithm;
class IDistanceCalculator;
class IVertexFit;
class IPVReFitter;

namespace LHCb {
class Particle;
class Vertex;
}  // namespace LHCb

class TupleToolApplyIsolation : public TupleToolBase,
                                virtual public IParticleTupleTool {
 public:
  TupleToolApplyIsolation( const std::string& type, const std::string& name,
                           const IInterface* parent );

  ~TupleToolApplyIsolation() override = default;

  StatusCode initialize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                   const std::string&, Tuples::Tuple& ) override;

 protected:
 private:
  Float_t opening, minipchi2, newfdchi2, oldfdchi2, ghostprob, trackchi2,
      deltafd, pt, ip, chi2, type, vertexchi2, Dst_PT, dummy;

  bool   isTrackInDecay( const LHCb::Track*, std::vector<const LHCb::Track*> );
  double getminipchi( const LHCb::Particle* );
  double getfdchi2( const LHCb::Track*, LHCb::Vertex );
  double getopening( const LHCb::Track*, const LHCb::Particle* );
  const LHCb::Vertex* originVertex( const LHCb::Particle*,
                                    const LHCb::Particle* ) const;

  IDVAlgorithm*            m_dva;
  IDistanceCalculator*     m_dist;
  const IVertexFit*        m_pVertexFit;
  IParticle2MCAssociator*  m_p2mcAssoc{};
  IPVReFitter*             m_pvReFitter{};
  TMVA::Reader*            m_Reader{};
  double                   m_deltaChi2;
  double                   m_Chi2;
  std::string              m_typeVertexFit;
  std::string              m_outputSuffix;
  std::string              m_weightsName;
  std::vector<std::string> m_inputParticles;
};

#endif  // TupleToolApplyIsolation_H
