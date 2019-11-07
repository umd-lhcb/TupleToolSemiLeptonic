#ifndef TupleToolDocas_H
#define TupleToolDocas_H

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h"
#include "LoKi/ChildSelector.h"

class IDVAlgorithm;
class IDistanceCalculator;

namespace LHCb {
class Particle;
}

/** @class TupleToolDocas TupleToolDocas.h
 *
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
 *  @author Paul Seyfert
 *  @date   2016-11-08
 */
class TupleToolDocas : public TupleToolBase,
                       virtual public IParticleTupleTool {
 public:
  /// Standard constructor
  TupleToolDocas( const std::string& type, const std::string& name,
                  const IInterface* parent );

  virtual ~TupleToolDocas(){};  ///< Destructor

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                   const std::string&, Tuples::Tuple& );

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

/**
 *
 *  DecayTreeTupleTool to compute DOCA between two particles in the decay chain
 *  Copyright (C) 2016  Paul Seyfert
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */
