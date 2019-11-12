#ifndef TupleToolDocas_H
#define TupleToolDocas_H 1

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h"
#include "LoKi/ChildSelector.h"

class IDVAlgorithm;
class IDistanceCalculator;

namespace LHCb {
class Particle;
}

/**
 * @class TupleToolDocas TupleToolDocas.h
 * @brief compute DOCAs between particles in the decay chain DecayTreeTuple
 *
 * Thanks to Vanya for the pointer to LoKi::Child::Selectors.
 *
 * a small example:
 * @code
 *
 * Bu2JpsiKTuple.Decay    = '[B- -> ^(J/psi(1S)-> ^mu+ ^mu-) ^K- ]CC'
 * Bu2JpsiKTuple.addBranches( {
 *     "muon_p"   : '[B- ->  (J/psi(1S)-> ^mu+  mu-)  K- ]CC'
 *    ,"muon_m"  : '[B- ->  (J/psi(1S)->  mu+ ^mu-)  K- ]CC'
 *    ,"kaon_m"  : '[B- ->  (J/psi(1S)->  mu+  mu-) ^K- ]CC'
 *    ,"Jpsi"    : '[B- -> ^(J/psi(1S)->  mu+  mu-)  K- ]CC'
 *    ,"Bu"      : '[B- ->  (J/psi(1S)->  mu+  mu-)  K- ]CC'
 *    } )
 *
 * budocas = Bu2JpsiKTuple.Bu.addTupleTool('TupleToolDocas')
 * budocas.Name       = [ "Kmu_OS", "Kmu_SS", "mumu" ]
 * budocas.Location1  = [ "[B- ->  (J/psi(1S)-> mu+  mu-)  ^K-]CC",
 *                        "[B- ->  (J/psi(1S)-> mu+  mu-)  ^K-]CC",
 *                        "[B- ->  (J/psi(1S)-> mu+  ^mu-)  K-]CC"]
 * budocas.Location2  = [ "[B- ->  (J/psi(1S)-> ^mu+  mu-)  K-]CC",
 *                        "[B- ->  (J/psi(1S)-> mu+  ^mu-)  K-]CC",
 *                        "[B- ->  (J/psi(1S)-> ^mu+  mu-)  K-]CC"]
 *
 * @endcode
 * creates the branches: Bu_DOCA_Kmu_OS, Bu_DOCA_Kmu_SS, Bu_DOCA_mumu
 * and: Bu_DOCACHI2_Kmu_OS, Bu_DOCACHI2_Kmu_SS, Bu_DOCACHI2_mumu
 *
 *
 * @sa DecayTreeTuple
 *
 * @author Paul Seyfert
 * @date   2016-11-08
 */
class TupleToolDocas : public TupleToolBase,
                       virtual public IParticleTupleTool {
 public:
  TupleToolDocas( const std::string& type, const std::string& name,
                  const IInterface* parent );

  ~TupleToolDocas() override = default;

  StatusCode initialize() override;
  StatusCode finalize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                   const std::string&, Tuples::Tuple& ) override;

 private:
  IDVAlgorithm*              m_dva;
  const IDistanceCalculator* m_dist;

  /// Well, I wanted a std::map<std::string,std::pair<std::string,std::string>>
  /// here, but Gaudi cannot parse that. Now I'm vulnerable against name
  /// collisions in m_name and lenght mismatches between the three vectors

  std::vector<std::string> m_locations1;
  std::vector<std::string> m_locations2;
  std::vector<std::string> m_name;

  std::map<std::string,
           std::pair<LoKi::Child::Selector*, LoKi::Child::Selector*>>
      m_childSelectors;
};

#endif  // TupleToolDocas_H
