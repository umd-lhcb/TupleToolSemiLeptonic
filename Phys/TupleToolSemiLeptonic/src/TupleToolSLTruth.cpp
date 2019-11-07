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

StatusCode TupleToolSLTruth::fill( const LHCb::Particle*,
                                   const LHCb::Particle* P,
                                   const std::string&    head,
                                   Tuples::Tuple&        tuple ) {
  const std::string prefix = fullName( head );

  bool test = true;

  const LHCb::MCParticle* mcp = getMCParticle( P );
  // pointer is ready, prepare the values
  // const int mcPid = ( mcp ? mcp->particleID().pid() : 0 );
  const SmartRefVector<LHCb::MCVertex>&    EndVertices = mcp->endVertices();
  SmartRefVector<MCVertex>::const_iterator itV         = EndVertices.begin();
  if ( EndVertices.size() > 1 ) ++itV;
  const LHCb::MCVertex*                  endV      = ( *itV );
  const SmartRefVector<LHCb::MCParticle> daughters = endV->products();
  // fill the tuple:
  LHCb::Particle::ConstVector tmp = P->daughtersVector();
  // ParticleVector tmp = m_particleDescendants->descendants(P);
  Gaudi::LorentzVector neutrinovector( 0, 0, 0, 0 );
  Gaudi::LorentzVector tauvector( 0, 0, 0, 0 );
  Gaudi::LorentzVector taunutauvector( 0, 0, 0, 0 );
  Gaudi::LorentzVector taunumuvector( 0, 0, 0, 0 );
  Gaudi::LorentzVector muonvector( 0, 0, 0, 0 );
  // Gaudi::LorentzVector hadronvector(0,0,0,0);
  std::vector<Gaudi::LorentzVector> hadronvectors;
  std::vector<int>                  hadronids;
  Gaudi::LorentzVector              mothervector( 0, 0, 0, 0 );
  // int hadronPID = 0;
  SmartRefVector<LHCb::MCParticle>::const_iterator itD;
  for ( itD = daughters.begin(); daughters.end() != itD; ++itD ) {
    if ( ( *itD )->particleID().abspid() == 22 ) continue;

    if ( abs( ( *itD )->particleID().pid() ) == 13 ) {
      muonvector = ( *itD )->momentum();
    }
    if ( abs( ( *itD )->particleID().pid() ) == 15 ) {
      tauvector = ( *itD )->momentum();
      const SmartRefVector<LHCb::MCVertex>& TauEndVertices =
          ( *itD )->endVertices();
      SmartRefVector<MCVertex>::const_iterator tauitV = TauEndVertices.begin();
      if ( TauEndVertices.size() > 1 ) ++tauitV;
      const LHCb::MCVertex*                  tauendV = ( *tauitV );
      const SmartRefVector<LHCb::MCParticle> taudaughters =
          tauendV->products();
      SmartRefVector<LHCb::MCParticle>::const_iterator itDtau;
      for ( itDtau = taudaughters.begin(); taudaughters.end() != itDtau;
            ++itDtau ) {
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
      hadronvectors.push_back( ( *itD )->momentum() );
      hadronids.push_back( ( *itD )->particleID().pid() );
    }
  }

  if ( mcp ) {
    mothervector = mcp->momentum();
  }

  Gaudi::LorentzVector munuvector = muonvector + neutrinovector;
  // Gaudi::LorentzVector sumvec = munuvector+charmvector;

  test &= tuple->column( prefix + "_True_Q2", munuvector.M2() );
  test &= tuple->column( prefix + "_TrueNeutrino_P", neutrinovector );
  test &= tuple->column( prefix + "_TrueTau_P", tauvector );
  test &= tuple->column( prefix + "_TrueTauNuMu_P", taunumuvector );
  test &= tuple->column( prefix + "_TrueTauNuTau_P", taunutauvector );
  for ( unsigned int i = 0; i < 3; i++ ) {
    std::stringstream sout;
    sout << i;
    if ( hadronvectors.size() > i ) {
      test &= tuple->column( prefix + "_TrueHadron_" + sout.str() + "_P",
                             hadronvectors[i] );
      test &= tuple->column( prefix + "_TrueHadron_" + sout.str() + "_ID",
                             hadronids[i] );
    } else {
      test &= tuple->column( prefix + "_TrueHadron_" + sout.str() + "_P",
                             Gaudi::LorentzVector( 0, 0, 0, 0 ) );
      test &= tuple->column( prefix + "_TrueHadron_" + sout.str() + "_ID", 0 );
    }
  }
  /*TLorentzVector Tcharmvector =
  TLorentzVector(charmvector.Px(),charmvector.Py(),charmvector.Pz(),charmvector.E());
  TLorentzVector Tmuonvector =
  TLorentzVector(muonvector.Px(),muonvector.Py(),muonvector.Pz(),muonvector.E());
  TLorentzVector Tneutrinovector =
  TLorentzVector(neutrinovector.Px(),neutrinovector.Py(),neutrinovector.Pz(),neutrinovector.E());
  //TLorentzVector mothervector = mcp->momentum();*/
  Gaudi::XYZVector     boostvector = munuvector.BoostToCM();
  Gaudi::LorentzVector restmother  = boostvec(
      mothervector, boostvector.X(), boostvector.Y(), boostvector.Z() );
  Gaudi::LorentzVector restmuon = boostvec( muonvector, boostvector.X(),
                                            boostvector.Y(), boostvector.Z() );
  double costhetal = -( restmother.Vect().Dot( restmuon.Vect() ) ) /
                     ( restmother.P() * restmuon.P() );
  test &= tuple->column( prefix + "_True_Costhetal", costhetal );
  // fill all requested MCTools
  for ( std::vector<IMCParticleTupleTool*>::const_iterator it =
            m_mcTools.begin();
        it != m_mcTools.end(); ++it ) {
    test &= ( *it )->fill( NULL, mcp, prefix, tuple );
  }

  return StatusCode( test );
}

//=============================================================================

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( TupleToolSLTruth )

//=============================================================================
