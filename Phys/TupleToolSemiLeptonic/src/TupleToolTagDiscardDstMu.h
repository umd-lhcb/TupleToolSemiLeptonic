#ifndef TupleToolTagDiscardDstMu_H
#define TupleToolTagDiscardDstMu_H 1

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h"  // Interface

class IDVAlgorithm;
class IDistanceCalculator;
class IVertexFit;
class ILifetimeFitter;

namespace LHCb {
class Particle;
class Vertex;
}  // namespace LHCb

class TupleToolTagDiscardDstMu : public TupleToolBase,
                                 virtual public IParticleTupleTool {
 public:
  TupleToolTagDiscardDstMu( const std::string& type, const std::string& name,
                            const IInterface* parent );

  ~TupleToolTagDiscardDstMu() override = default;  ///< Destructor

  StatusCode initialize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                   const std::string&, Tuples::Tuple& ) override;

 private:
  const LHCb::Vertex* originVertex( const LHCb::Particle*,
                                    const LHCb::Particle* ) const;

 private:
  IDVAlgorithm*              m_dva;
  const IDistanceCalculator* m_dist;
  const IVertexFit*          m_pVertexFit;
  ILifetimeFitter*           m_ltfit;
  std::string                m_outputSuffix;
  int                        m_pri1;
  int                        m_pri2;
  int                        m_pri3;
  int                        m_sec1;
  int                        m_sec2;
  int                        m_sec3;
  int                        m_ter1;
  int                        m_ter2;
  int                        m_ter3;
  int                        m_priDOCA;
  int                        m_secDOCA;
  int                        m_terDOCA;
  int                        m_priexDOCA1;
  int                        m_secexDOCA1;
  int                        m_terexDOCA1;
  int                        m_priexDOCA2;
  int                        m_secexDOCA2;
  int                        m_terexDOCA2;
  int                        m_TimePID;
  bool                       m_Dstar;
  bool                       m_discSlow;
  std::string                m_typeVertexFit;
  std::vector<std::string>   m_inputParticles;
};

#endif  // TupleToolTagDiscardDstMu_H
