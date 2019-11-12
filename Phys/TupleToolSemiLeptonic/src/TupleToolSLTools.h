#ifndef TupleToolSLTools_H
#define TupleToolSLTools_H 1

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h"  // Interface
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TrackInterfaces/ITrackVertexer.h"

class IDVAlgorithm;
class IDistanceCalculator;
class IVertexFit;
class ILifetimeFitter;
class IPVReFitter;

namespace LHCb {
class Particle;
class Vertex;
}  // namespace LHCb

class TupleToolSLTools : public TupleToolBase,
                         virtual public IParticleTupleTool {
 public:
  TupleToolSLTools( const std::string& type, const std::string& name,
                    const IInterface* parent );

  ~TupleToolSLTools() override = default;

  StatusCode initialize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                   const std::string&, Tuples::Tuple& ) override;

 private:
  double Mcorr( TLorentzVector, TVector3 );

  std::vector<double> Mcorr_errors( TVector3 v1, TVector3 v2, TLorentzVector p,
                                    Gaudi::SymMatrix7x7 cov1,
                                    Gaudi::SymMatrix3x3 cov2 );
  std::vector<TLorentzVector> recoNu( TLorentzVector Y, TVector3 M_dirn,
                                      double mass );
  double dpnuhighdpperp( double pperp, TLorentzVector visible );
  double dpnulowdpperp( double pperp, TLorentzVector visible );
  double dq2dpnu( double, TLorentzVector, TVector3, TLorentzVector );
  std::vector<double> q2_errors( TVector3, TVector3, TLorentzVector,
                                 TLorentzVector, Gaudi::SymMatrix7x7,
                                 Gaudi::SymMatrix3x3 );
  StatusCode          fillVertex( const LHCb::VertexBase* vtx,
                                  const std::string&      vtx_name,
                                  Tuples::Tuple&          tuple ) const;

  IDVAlgorithm*              m_dva;
  const IDistanceCalculator* m_dist;
  const IVertexFit*          m_pVertexFit;
  std::string                m_typeVertexFit;
  double                     m_Bmass;
  bool                       m_vcov;
  bool                       m_momcov;
};

#endif  // TupleToolSLTools_H
