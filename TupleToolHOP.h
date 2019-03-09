#ifndef TUPLETOOLHOP_H
#define TUPLETOOLHOP_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/IParticleTupleTool.h"  // Interface
#include "Kernel/ParticleProperty.h"

class IDVAlgorithm;

namespace LHCb {
class Particle;
class VertexBase;
}  // namespace LHCb

/** @class TupleToolHOP TupleToolHOP.h
 *
 * \brief Fill geometry related information for DecayTreeTuple
 *
 * - mother_HOP : HOP correction, the ratio between the Pt (with respecto to the
 * mothers direction of flight) of electrons and non-electrons
 * - mother_HOP_MASS : Mother mass obtained by scalling the electrons P by the
 * hop factor
 * - mother_ELECTRON_MASS : Mass obtained by summing the four-momenta of all
 * electrons
 *
 *
 *
 * \sa DecayTreeTuple
 *
 *  @author Vinicius Franco, Carla Marin Benito
 *  @date   2016-29-06
 */

struct ClassifyParticles_Tuple {
  int isElectron;       // 1 == electron, 0 == other, 2 == no type
  int hasSameChildren;  // 1 == yes, 0 == no, 2 == no children
};

class TupleToolHOP : public TupleToolBase, virtual public IParticleTupleTool {
 public:
  /// Standard constructor
  TupleToolHOP(const std::string& type, const std::string& name,
               const IInterface* parent);

  virtual ~TupleToolHOP(){};  ///< Destructor

  StatusCode initialize() override;

  StatusCode fill(const LHCb::Particle*, const LHCb::Particle*,
                  const std::string&, Tuples::Tuple&) override;

 private:
  float HOPProjectMomentum(const LHCb::Particle* top,
                           const Gaudi::LorentzVector* part_four_mom) const;

  ClassifyParticles_Tuple* ClassifyParticles(
      const LHCb::Particle& top,
      SmartRefVector<LHCb::Particle>& electronContainer,
      SmartRefVector<LHCb::Particle>& nonElectronContainer) const;

  ClassifyParticles_Tuple* ClassifyParticles_Merge(
      const SmartRefVector<LHCb::Particle>& dau,
      //      SmartRefVector<ClassifyParticles_Tuple>& t_list,
      std::vector<ClassifyParticles_Tuple*>& t_list, bool are_equal,
      SmartRefVector<LHCb::Particle>& electronContainer,
      SmartRefVector<LHCb::Particle>& nonElectronContainer) const;

 private:
  // bool m_fillMother;
  IDVAlgorithm* m_dva;
  LHCb::IParticlePropertySvc* m_ppSvc;
  float m_electronMassSquared = 0;
};
#endif  // TUPLETOOLHOP_H
