// $Id: TupleToolSLTruth.cpp,v 1.17 2010-01-26 15:39:26 rlambert Exp $
// Include files
#include "gsl/gsl_sys.h"
// from Gaudi
#include "GaudiKernel/PhysicalConstants.h"
// local
#include "TupleToolSLTruth.h"

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/MCParticle.h"
#include "Event/Particle.h"
#include "TLorentzVector.h"
// kernel
#include "Kernel/IParticle2MCAssociator.h"
//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolSLTruth
//
// 2007-11-07 : Jeremie Borel
//-----------------------------------------------------------------------------

using namespace LHCb;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolSLTruth::TupleToolSLTruth( const std::string& type,
                                    const std::string& name,
                                    const IInterface*  parent )
    : TupleToolBase( type, name, parent ),
      m_toolList( 1, "MCTupleToolKinematic" ) {
  // interface
  declareInterface<IParticleTupleTool>( this );
  // The names of MCTupleTools to use on the associated mcp
  declareProperty( "ToolList", m_toolList );
  // MC associators to try, in order
  m_p2mcAssocTypes.push_back( "DaVinciSmartAssociator" );
  m_p2mcAssocTypes.push_back( "MCMatchObjP2MCRelator" );
  declareProperty( "IP2MCPAssociatorTypes", m_p2mcAssocTypes );
}

//=============================================================================

StatusCode TupleToolSLTruth::initialize() {
  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;
  m_particleDescendants =
      tool<IParticleDescendants>( "ParticleDescendants", this );
  // the MC associators
  m_p2mcAssocs.clear();
  for ( std::vector<std::string>::const_iterator iMCAss =
            m_p2mcAssocTypes.begin();
        iMCAss != m_p2mcAssocTypes.end(); ++iMCAss ) {
    m_p2mcAssocs.push_back( tool<IParticle2MCAssociator>( *iMCAss, this ) );
  }
  if ( m_p2mcAssocs.empty() ) {
    return Error( "No MC associators configured" );
  }

  // remove duplicate tools from the list
  std::sort( m_toolList.begin(), m_toolList.end() );
  m_toolList.erase( std::unique( m_toolList.begin(), m_toolList.end() ),
                    m_toolList.end() );

  // initialise the tuple tools
  for ( std::vector<std::string>::const_iterator it = m_toolList.begin();
        m_toolList.end() != it; ++it ) {
    if ( msgLevel( MSG::VERBOSE ) )
      verbose() << "Adding the tool " << *it << endmsg;
    IMCParticleTupleTool* aTool = tool<IMCParticleTupleTool>( *it, this );
    if ( aTool ) {
      m_mcTools.push_back( aTool );
    } else {
      Warning( "There was a problem retrieving " + *it +
               " , this tool will be ignored" )
          .ignore();
    }
  }

  if ( msgLevel( MSG::VERBOSE ) ) {
    verbose() << "Completed TupleTool intialisation, " << m_mcTools.size()
              << " tools added " << endmsg;
  }

  return sc;
}

//=============================================================================

bool TupleToolSLTruth::isCharmHadron( int PID ) {
  bool isCharm = false;
  if ( PID < 100 ) return false;
  while ( PID ) {
    if ( abs( PID ) % 10 == 4 ) {
      isCharm = true;
    }
    PID /= 10;
  }
  return isCharm;
}

bool TupleToolSLTruth::isBeautyHadron( int PID ) {
  bool isBeauty = false;
  while ( PID ) {
    if ( abs( PID ) % 10 == 5 ) {
      isBeauty = true;
    }
    PID /= 10;
  }
  return isBeauty;
}

const LHCb::MCParticle* TupleToolSLTruth::getMCParticle(
    const LHCb::Particle* P ) {
  const LHCb::MCParticle* mcp( NULL );
  if ( P ) {
    // assignedPid = P->particleID().pid();
    if ( msgLevel( MSG::VERBOSE ) )
      verbose() << "Getting related MCP to " << P << endmsg;
    for ( std::vector<IParticle2MCAssociator*>::const_iterator iMCAss =
              m_p2mcAssocs.begin();
          iMCAss != m_p2mcAssocs.end(); ++iMCAss ) {
      mcp = ( *iMCAss )->relatedMCP( P );
      if ( mcp ) break;
    }
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Got mcp " << mcp << endmsg;
  }
  return mcp;
}

bool TupleToolSLTruth::motherIsBeauty( LHCb::MCParticle*       MCPart,
                                       const LHCb::MCParticle* mcp ) {
  bool motherBeauty = false;
  if ( mcp && MCPart && MCPart->mother() ) {
    if ( MCPart->mother() == mcp ) {
      motherBeauty = true;
    }
  }
  return motherBeauty;
}

Gaudi::LorentzVector boostvec( Gaudi::LorentzVector vec, Double_t bx,
                               Double_t by, Double_t bz ) {
  Gaudi::LorentzVector newvec;
  // Boost this Lorentz vector
  Double_t b2     = bx * bx + by * by + bz * bz;
  Double_t gamma  = 1.0 / TMath::Sqrt( 1.0 - b2 );
  Double_t bp     = bx * vec.X() + by * vec.Y() + bz * vec.Z();
  Double_t gamma2 = b2 > 0 ? ( gamma - 1.0 ) / b2 : 0.0;
  newvec.SetPx( vec.X() + gamma2 * bp * bx + gamma * bx * vec.T() );
  newvec.SetPy( vec.Y() + gamma2 * bp * by + gamma * by * vec.T() );
  newvec.SetPz( vec.Z() + gamma2 * bp * bz + gamma * bz * vec.T() );
  newvec.SetE( gamma * ( vec.T() + bp ) );

  return newvec;
}

StatusCode TupleToolSLTruth::getMCTruePi0s( const LHCb::MCParticle* MCP,
                                            std::vector<double>* MCTruePi0PX_,
                                            std::vector<double>* MCTruePi0PY_,
                                            std::vector<double>* MCTruePi0PZ_,
                                            std::vector<double>* MCTruePi0E_,
                                            std::vector<int>* MCTruePi0MID_ ) {
  if ( MCP ) {
    if ( abs( MCP->particleID().pid() ) == 111 ) {
      Gaudi::LorentzVector MCPmomentum = MCP->momentum();
      if ( MCPmomentum.E() > 0 && MCP->mother() ) {
        MCTruePi0PX_->push_back( MCPmomentum.X() );
        MCTruePi0PY_->push_back( MCPmomentum.Y() );
        MCTruePi0PZ_->push_back( MCPmomentum.Z() );
        MCTruePi0E_->push_back( MCPmomentum.E() );
        MCTruePi0MID_->push_back( MCP->mother()->particleID().pid() );
      }
    } else {
      if ( !MCP->endVertices().empty() &&
           abs( MCP->particleID().pid() ) > 100 &&
           abs( MCP->particleID().pid() ) != 211 and
           abs( MCP->particleID().pid() ) != 321 ) {
        const SmartRefVector<LHCb::MCVertex>& MCPEndVertices =
            MCP->endVertices();
        SmartRefVector<MCVertex>::const_iterator MCPitV =
            MCPEndVertices.begin();
        if ( MCPEndVertices.size() > 1 ) ++MCPitV;
        const LHCb::MCVertex* MCPendV = ( *MCPitV );
        if ( MCPendV ) {
          const SmartRefVector<LHCb::MCParticle> MCPdaughters =
              MCPendV->products();
          if ( !MCPdaughters.empty() ) {
            SmartRefVector<LHCb::MCParticle>::const_iterator itDMCP;
            for ( itDMCP = MCPdaughters.begin(); MCPdaughters.end() != itDMCP;
                  ++itDMCP ) {
              if ( ( *itDMCP ) ) {
                getMCTruePi0s( ( *itDMCP ), MCTruePi0PX_, MCTruePi0PY_,
                               MCTruePi0PZ_, MCTruePi0E_, MCTruePi0MID_ );
              }
            }
          }
        }
      }
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode TupleToolSLTruth::getMCTrueGammas(
    const LHCb::MCParticle* MCP, std::vector<double>* MCTrueGammaPX_,
    std::vector<double>* MCTrueGammaPY_, std::vector<double>* MCTrueGammaPZ_,
    std::vector<double>* MCTrueGammaE_, std::vector<int>* MCTrueGammaMID_,
    std::vector<int>* MCTrueGammaIsRad_ ) {
  if ( MCP ) {
    int MCPid = abs( MCP->particleID().pid() );
    if ( MCP->mother() ) {
      if ( MCPid == 22 ) {
        Gaudi::LorentzVector MCPmomentum = MCP->momentum();
        if ( MCPmomentum.E() > 0 ) {
          MCTrueGammaPX_->push_back( MCPmomentum.X() );
          MCTrueGammaPY_->push_back( MCPmomentum.Y() );
          MCTrueGammaPZ_->push_back( MCPmomentum.Z() );
          MCTrueGammaE_->push_back( MCPmomentum.E() );
          MCTrueGammaMID_->push_back( MCP->mother()->particleID().pid() );
          MCTrueGammaIsRad_->push_back(
              -1 );  // I still don't know how to fill this attribute properly,
                     // skip it
        }
      } else {
        if ( !MCP->endVertices().empty() ) {
          const SmartRefVector<LHCb::MCVertex>& MCPEndVertices =
              MCP->endVertices();
          SmartRefVector<MCVertex>::const_iterator MCPitV =
              MCPEndVertices.begin();
          if ( MCPEndVertices.size() > 1 ) ++MCPitV;
          const LHCb::MCVertex* MCPendV = ( *MCPitV );
          if ( MCPendV ) {
            const SmartRefVector<LHCb::MCParticle> MCPdaughters =
                MCPendV->products();
            if ( !MCPdaughters.empty() ) {
              SmartRefVector<LHCb::MCParticle>::const_iterator itDMCP;
              bool inspectdaughters = 1;
              if ( MCPid == 211 || MCPid == 321 || MCPid == 13 ) {
                inspectdaughters = 0;
              }
              if ( inspectdaughters == 1 ) {
                for ( itDMCP = MCPdaughters.begin();
                      MCPdaughters.end() != itDMCP; ++itDMCP ) {
                  if ( ( *itDMCP ) ) {
                    getMCTrueGammas( ( *itDMCP ), MCTrueGammaPX_,
                                     MCTrueGammaPY_, MCTrueGammaPZ_,
                                     MCTrueGammaE_, MCTrueGammaMID_,
                                     MCTrueGammaIsRad_ );
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode TupleToolSLTruth::fill( const LHCb::Particle*,
                                   const LHCb::Particle* P,
                                   const std::string&    head,
                                   Tuples::Tuple&        tuple ) {
  const std::string prefix = fullName( head );

  bool test = true;

  Gaudi::LorentzVector neutrinovector( 0, 0, 0, 0 );
  Gaudi::LorentzVector tauvector( 0, 0, 0, 0 );
  Gaudi::LorentzVector taunutauvector( 0, 0, 0, 0 );
  Gaudi::LorentzVector taunumuvector( 0, 0, 0, 0 );
  Gaudi::LorentzVector muonvector( 0, 0, 0, 0 );
  Gaudi::LorentzVector munuvector( 0, 0, 0, 0 );
  // Gaudi::LorentzVector hadronvector(0,0,0,0);
  std::vector<Gaudi::LorentzVector> hadronvectors_D;
  std::vector<Gaudi::LorentzVector> hadronvectors_D1GD;
  std::vector<Gaudi::LorentzVector> hadronvectors_D2GD;
  std::vector<Gaudi::LorentzVector> hadronvectors_D3GD;
  std::vector<Gaudi::LorentzVector> hadronvectors_D4GD;
  std::vector<int>                  hadronids_D;
  std::vector<int>                  hadronids_D1GD;
  std::vector<int>                  hadronids_D2GD;
  std::vector<int>                  hadronids_D3GD;
  std::vector<int>                  hadronids_D4GD;
  Gaudi::LorentzVector              mothervector( 0, 0, 0, 0 );
  double                            costhetal = 10.;
  bool                              is_tau    = false;

  int                 MCTruePi0Multiplicity   = 0;
  int                 MCTrueGammaMultiplicity = 0;
  std::vector<double> MCTruePi0PX, MCTruePi0PY, MCTruePi0PZ, MCTruePi0E,
      MCTrueGammaPX, MCTrueGammaPY, MCTrueGammaPZ, MCTrueGammaE;
  std::vector<int> MCTruePi0MID, MCTrueGammaMID, MCTrueGammaIsRad;
  double           MCTrueCombGammaPX = 0.;
  double           MCTrueCombGammaPY = 0.;
  double           MCTrueCombGammaPZ = 0.;
  double           MCTrueCombGammaE  = 0.;
  double           MCTrueCombGammaPT = 0.;

  const LHCb::MCParticle* mcp = getMCParticle( P );
  if ( mcp ) {
    if ( !mcp->endVertices().empty() ) {
      // if (abs(mcp->particleID().pid())==511 ||
      // abs(mcp->particleID().pid())==521 || abs(mcp->particleID().pid())==531
      // || abs(mcp->particleID().pid())==5221) {
      if ( isBeautyHadron( abs( mcp->particleID().pid() ) ) ) {
        getMCTruePi0s( mcp, &MCTruePi0PX, &MCTruePi0PY, &MCTruePi0PZ,
                       &MCTruePi0E, &MCTruePi0MID );
        getMCTrueGammas( mcp, &MCTrueGammaPX, &MCTrueGammaPY, &MCTrueGammaPZ,
                         &MCTrueGammaE, &MCTrueGammaMID, &MCTrueGammaIsRad );
      }

      // pointer is ready, prepare the values
      // const int mcPid = ( mcp ? mcp->particleID().pid() : 0 );
      const SmartRefVector<LHCb::MCVertex>& EndVertices = mcp->endVertices();
      SmartRefVector<MCVertex>::const_iterator itV      = EndVertices.begin();
      SmartRefVector<MCVertex>::const_iterator itVD;
      if ( EndVertices.size() > 1 ) ++itV;
      const LHCb::MCVertex*                  endV      = ( *itV );
      const SmartRefVector<LHCb::MCParticle> daughters = endV->products();
      SmartRefVector<LHCb::MCParticle>       granddaughters;
      // fill the tuple:
      LHCb::Particle::ConstVector tmp = P->daughtersVector();
      // ParticleVector tmp = m_particleDescendants->descendants(P);
      // int hadronPID = 0;
      SmartRefVector<LHCb::MCParticle>::const_iterator itD;
      int                                              D_counter = 0;
      SmartRefVector<LHCb::MCParticle>::const_iterator itGD;
      for ( itD = daughters.begin(); daughters.end() != itD; ++itD ) {
        if ( ( *itD )->particleID().abspid() == 22 ) continue;

        if ( abs( ( *itD )->particleID().pid() ) == 13 ) {
          muonvector = ( *itD )->momentum();
        }
        if ( ( abs( ( *itD )->particleID().pid() ) == 15 ) &&
             !( *itD )->endVertices().empty() ) {
          is_tau    = true;
          tauvector = ( *itD )->momentum();
          const SmartRefVector<LHCb::MCVertex>& TauEndVertices =
              ( *itD )->endVertices();
          SmartRefVector<MCVertex>::const_iterator tauitV =
              TauEndVertices.begin();
          if ( TauEndVertices.size() > 1 ) ++tauitV;
          const LHCb::MCVertex*                  tauendV = ( *tauitV );
          const SmartRefVector<LHCb::MCParticle> taudaughters =
              tauendV->products();
          SmartRefVector<LHCb::MCParticle>::const_iterator itDtau;
          for ( itDtau = taudaughters.begin(); taudaughters.end() != itDtau;
                ++itDtau ) {
            if ( abs( ( *itDtau )->particleID().pid() ) == 13 ) {
              muonvector = ( *itDtau )->momentum();
            }
            if ( abs( ( *itDtau )->particleID().pid() ) == 14 ) {
              taunumuvector = ( *itDtau )->momentum();
            }
            if ( abs( ( *itDtau )->particleID().pid() ) == 16 ) {
              taunutauvector = ( *itDtau )->momentum();
            }
          }
        }
        if ( ( abs( ( *itD )->particleID().pid() ) == 14 ) ||
             ( abs( ( *itD )->particleID().pid() ) == 16 ) ) {
          neutrinovector = ( *itD )->momentum();
        }
        if ( abs( ( *itD )->particleID().pid() ) > 100 ) {
          hadronvectors_D.push_back( ( *itD )->momentum() );
          hadronids_D.push_back( ( *itD )->particleID().pid() );

          if ( !( *itD )->endVertices().empty() ) {
            itVD = ( *itD )->endVertices().begin();
            if ( ( *itD )->endVertices().size() > 1 ) ++itVD;
            granddaughters = ( *itVD )->products();

            for ( itGD = granddaughters.begin(); granddaughters.end() != itGD;
                  ++itGD ) {
              if ( D_counter == 0 ) {
                hadronvectors_D1GD.push_back( ( *itGD )->momentum() );
                hadronids_D1GD.push_back( ( *itGD )->particleID().pid() );
              } else if ( D_counter == 1 ) {
                hadronvectors_D2GD.push_back( ( *itGD )->momentum() );
                hadronids_D2GD.push_back( ( *itGD )->particleID().pid() );
              } else if ( D_counter == 2 ) {
                hadronvectors_D3GD.push_back( ( *itGD )->momentum() );
                hadronids_D3GD.push_back( ( *itGD )->particleID().pid() );
              } else if ( D_counter == 3 ) {
                hadronvectors_D4GD.push_back( ( *itGD )->momentum() );
                hadronids_D4GD.push_back( ( *itGD )->particleID().pid() );
              }
            }
            D_counter += 1;
          }
        }
      }

      mothervector = mcp->momentum();
      if ( is_tau )
        munuvector = tauvector + neutrinovector;
      else
        munuvector = muonvector + neutrinovector;
      // Gaudi::LorentzVector sumvec = munuvector+charmvector;
      Gaudi::XYZVector     boostvector = munuvector.BoostToCM();
      Gaudi::LorentzVector restmother  = boostvec(
          mothervector, boostvector.X(), boostvector.Y(), boostvector.Z() );
      Gaudi::LorentzVector restmuon = boostvec(
          muonvector, boostvector.X(), boostvector.Y(), boostvector.Z() );
      if ( restmother.P() * restmuon.P() > 0 ) {
        costhetal = -( restmother.Vect().Dot( restmuon.Vect() ) ) /
                    ( restmother.P() * restmuon.P() );
      }
    }
  }

  test &= tuple->column( prefix + "_True_IsTauDecay", is_tau );
  test &= tuple->column( prefix + "_True_Q2", munuvector.M2() );
  test &= tuple->column( prefix + "_TrueNeutrino_P", neutrinovector );
  test &= tuple->column( prefix + "_TrueMu_P", muonvector );
  test &= tuple->column( prefix + "_TrueTau_P", tauvector );
  test &= tuple->column( prefix + "_TrueTauNuMu_P", taunumuvector );
  test &= tuple->column( prefix + "_TrueTauNuTau_P", taunutauvector );

  for ( unsigned int i = 0; i < 4; i++ ) {
    std::stringstream sout;
    sout << i;
    if ( hadronvectors_D.size() > i ) {
      test &= tuple->column( prefix + "_TrueHadron_D" + sout.str() + "_P",
                             hadronvectors_D[i] );
      test &= tuple->column( prefix + "_TrueHadron_D" + sout.str() + "_ID",
                             hadronids_D[i] );
    } else {
      test &= tuple->column( prefix + "_TrueHadron_D" + sout.str() + "_P",
                             Gaudi::LorentzVector( 0, 0, 0, 0 ) );
      test &=
          tuple->column( prefix + "_TrueHadron_D" + sout.str() + "_ID", 0 );
    }
  }
  for ( unsigned int j = 0; j < 4; j++ ) {
    std::stringstream sout1;
    sout1 << j;
    if ( hadronvectors_D1GD.size() > j ) {
      test &= tuple->column( prefix + "_TrueHadron_D0_GD" + sout1.str() + "_P",
                             hadronvectors_D1GD[j] );
      test &=
          tuple->column( prefix + "_TrueHadron_D0_GD" + sout1.str() + "_ID",
                         hadronids_D1GD[j] );
    } else {
      test &= tuple->column( prefix + "_TrueHadron_D0_GD" + sout1.str() + "_P",
                             Gaudi::LorentzVector( 0, 0, 0, 0 ) );
      test &= tuple->column(
          prefix + "_TrueHadron_D0_GD" + sout1.str() + "_ID", 0 );
    }
  }
  for ( unsigned int j = 0; j < 4; j++ ) {
    std::stringstream sout2;
    sout2 << j;
    if ( hadronvectors_D2GD.size() > j ) {
      test &= tuple->column( prefix + "_TrueHadron_D1_GD" + sout2.str() + "_P",
                             hadronvectors_D2GD[j] );
      test &=
          tuple->column( prefix + "_TrueHadron_D1_GD" + sout2.str() + "_ID",
                         hadronids_D2GD[j] );
    } else {
      test &= tuple->column( prefix + "_TrueHadron_D1_GD" + sout2.str() + "_P",
                             Gaudi::LorentzVector( 0, 0, 0, 0 ) );
      test &= tuple->column(
          prefix + "_TrueHadron_D1_GD" + sout2.str() + "_ID", 0 );
    }
  }
  for ( unsigned int j = 0; j < 4; j++ ) {
    std::stringstream sout3;
    sout3 << j;
    if ( hadronvectors_D3GD.size() > j ) {
      test &= tuple->column( prefix + "_TrueHadron_D2_GD" + sout3.str() + "_P",
                             hadronvectors_D3GD[j] );
      test &=
          tuple->column( prefix + "_TrueHadron_D2_GD" + sout3.str() + "_ID",
                         hadronids_D3GD[j] );
    } else {
      test &= tuple->column( prefix + "_TrueHadron_D2_GD" + sout3.str() + "_P",
                             Gaudi::LorentzVector( 0, 0, 0, 0 ) );
      test &= tuple->column(
          prefix + "_TrueHadron_D2_GD" + sout3.str() + "_ID", 0 );
    }
  }
  for ( unsigned int j = 0; j < 4; j++ ) {
    std::stringstream sout4;
    sout4 << j;
    if ( hadronvectors_D4GD.size() > j ) {
      test &= tuple->column( prefix + "_TrueHadron_D3_GD" + sout4.str() + "_P",
                             hadronvectors_D4GD[j] );
      test &=
          tuple->column( prefix + "_TrueHadron_D3_GD" + sout4.str() + "_ID",
                         hadronids_D4GD[j] );
    } else {
      test &= tuple->column( prefix + "_TrueHadron_D3_GD" + sout4.str() + "_P",
                             Gaudi::LorentzVector( 0, 0, 0, 0 ) );
      test &= tuple->column(
          prefix + "_TrueHadron_D3_GD" + sout4.str() + "_ID", 0 );
    }
  }

  test &= tuple->column( prefix + "_True_Costhetal", costhetal );

  MCTruePi0Multiplicity   = MCTruePi0E.size();
  MCTrueGammaMultiplicity = MCTrueGammaE.size();
  std::vector<double> MCTruePi0PX_copy, MCTruePi0PY_copy, MCTruePi0PZ_copy,
      MCTruePi0E_copy, MCTrueGammaPX_copy, MCTrueGammaPY_copy,
      MCTrueGammaPZ_copy, MCTrueGammaE_copy;
  std::vector<int> MCTruePi0MID_copy, MCTrueGammaMID_copy,
      MCTrueGammaIsRad_copy;
  if ( MCTruePi0Multiplicity != 0 ) {
    for ( int i = 0; i < MCTruePi0Multiplicity; i++ ) {
      MCTruePi0PX_copy.push_back( MCTruePi0PX[i] );
      MCTruePi0PY_copy.push_back( MCTruePi0PY[i] );
      MCTruePi0PZ_copy.push_back( MCTruePi0PZ[i] );
      MCTruePi0E_copy.push_back( MCTruePi0E[i] );
      MCTruePi0MID_copy.push_back( MCTruePi0MID[i] );
    }
  } else {
    MCTruePi0PX_copy  = { 0 };
    MCTruePi0PY_copy  = { 0 };
    MCTruePi0PZ_copy  = { 0 };
    MCTruePi0E_copy   = { 0 };
    MCTruePi0MID_copy = { 0 };
  }
  if ( MCTrueGammaMultiplicity != 0 ) {
    for ( int i = 0; i < MCTrueGammaMultiplicity; i++ ) {
      MCTrueGammaPX_copy.push_back( MCTrueGammaPX[i] );
      MCTrueGammaPY_copy.push_back( MCTrueGammaPY[i] );
      MCTrueGammaPZ_copy.push_back( MCTrueGammaPZ[i] );
      MCTrueGammaE_copy.push_back( MCTrueGammaE[i] );
      MCTrueGammaMID_copy.push_back( MCTrueGammaMID[i] );
      MCTrueGammaIsRad_copy.push_back( MCTrueGammaIsRad[i] );
      MCTrueCombGammaPX += MCTrueGammaPX[i];
      MCTrueCombGammaPY += MCTrueGammaPY[i];
      MCTrueCombGammaPZ += MCTrueGammaPZ[i];
      MCTrueCombGammaE += MCTrueGammaE[i];
    }
  } else {
    MCTrueGammaPX_copy    = { 0 };
    MCTrueGammaPY_copy    = { 0 };
    MCTrueGammaPZ_copy    = { 0 };
    MCTrueGammaE_copy     = { 0 };
    MCTrueGammaMID_copy   = { 0 };
    MCTrueGammaIsRad_copy = { 0 };
  }
  MCTrueCombGammaPT = TMath::Sqrt( MCTrueCombGammaPX * MCTrueCombGammaPX +
                                   MCTrueCombGammaPY * MCTrueCombGammaPY );

  const int max_neutral_candidates = 1000;

  test &= tuple->column( prefix + "_MCTrue_pi0_mult", MCTruePi0Multiplicity );
  test &= tuple->farray( prefix + "_MCTrue_pi0_PX", MCTruePi0PX_copy,
                         prefix + "_MCTrue_pi0_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray( prefix + "_MCTrue_pi0_PY", MCTruePi0PY_copy,
                         prefix + "_MCTrue_pi0_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray( prefix + "_MCTrue_pi0_PZ", MCTruePi0PZ_copy,
                         prefix + "_MCTrue_pi0_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray( prefix + "_MCTrue_pi0_E", MCTruePi0E_copy,
                         prefix + "_MCTrue_pi0_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray( prefix + "_MCTrue_pi0_mother_ID", MCTruePi0MID_copy,
                         prefix + "_MCTrue_pi0_ArrayLength",
                         max_neutral_candidates );
  test &=
      tuple->column( prefix + "_MCTrue_gamma_mult", MCTrueGammaMultiplicity );
  test &= tuple->farray( prefix + "_MCTrue_gamma_PX", MCTrueGammaPX_copy,
                         prefix + "_MCTrue_gamma_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray( prefix + "_MCTrue_gamma_PY", MCTrueGammaPY_copy,
                         prefix + "_MCTrue_gamma_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray( prefix + "_MCTrue_gamma_PZ", MCTrueGammaPZ_copy,
                         prefix + "_MCTrue_gamma_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray( prefix + "_MCTrue_gamma_E", MCTrueGammaE_copy,
                         prefix + "_MCTrue_gamma_ArrayLength",
                         max_neutral_candidates );
  test &= tuple->farray(
      prefix + "_MCTrue_gamma_mother_ID", MCTrueGammaMID_copy,
      prefix + "_MCTrue_gamma_ArrayLength", max_neutral_candidates );
  test &= tuple->farray(
      prefix + "_MCTrue_gamma_is_radiative", MCTrueGammaIsRad_copy,
      prefix + "_MCTrue_gamma_ArrayLength", max_neutral_candidates );
  test &= tuple->column( prefix + "_MCTrue_combgamma_PX", MCTrueCombGammaPX );
  test &= tuple->column( prefix + "_MCTrue_combgamma_PY", MCTrueCombGammaPY );
  test &= tuple->column( prefix + "_MCTrue_combgamma_PZ", MCTrueCombGammaPZ );
  test &= tuple->column( prefix + "_MCTrue_combgamma_PT", MCTrueCombGammaPT );
  test &= tuple->column( prefix + "_MCTrue_combgamma_E", MCTrueCombGammaE );

  if ( mcp ) {
    if ( !mcp->endVertices().empty() ) {
      // fill all requested MCTools
      for ( std::vector<IMCParticleTupleTool*>::const_iterator it =
                m_mcTools.begin();
            it != m_mcTools.end(); ++it ) {
        test &= ( *it )->fill( NULL, mcp, prefix, tuple );
      }
    }
  }

  return StatusCode( test );
}

//=============================================================================

DECLARE_COMPONENT( TupleToolSLTruth )

//=============================================================================
