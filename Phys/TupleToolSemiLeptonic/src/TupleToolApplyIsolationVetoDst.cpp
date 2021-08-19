// local
#include "TupleToolApplyIsolationVetoDst.h"

// from Gaudi
#include "GaudiAlg/Tuple.h"

// from Phys
#include "Kernel/GetIDVAlgorithm.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IDistanceCalculator.h"
#include "Kernel/IVertexFit.h"

// from LHCb
#include "Event/MCParticle.h"
#include "Event/Particle.h"
#include "Kernel/IPVReFitter.h"
#include "Linker/LinkerTable.h"
#include "TrackInterfaces/ITrackVertexer.h"

// from ROOT
#include "TFile.h"
#include "TPluginManager.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"

#include <functional>
#include <string>

DECLARE_COMPONENT( TupleToolApplyIsolationVetoDst )

using namespace LHCb;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolApplyIsolationVetoDst::TupleToolApplyIsolationVetoDst(
    const std::string& type, const std::string& name,
    const IInterface* parent )
    : TupleToolBase( type, name, parent ),
      m_dva( nullptr ),
      m_dist( nullptr ),
      m_pVertexFit( nullptr ) {
  declareInterface<IParticleTupleTool>( this );

  m_inputParticles.emplace_back( "/Event/Phys/StdAllNoPIDsPions" );
  m_inputParticles.emplace_back( "/Event/Phys/StdNoPIDsUpPions" );
  m_inputParticles.emplace_back( "Phys/StdNoPIDsVeloPions" );

  declareProperty( "VertexFit", m_typeVertexFit = "default" );
  declareProperty( "InputParticles", m_inputParticles );
  declareProperty( "OutputSuffix", m_outputSuffix = "" );
  declareProperty( "WeightsFile", m_weightsName = "weights.xml" );
  declareProperty( "NWrite", m_nWrite = 3 );
  declareProperty( "TrueIDs", m_trueID = false );
  declareProperty( "Verbose", m_verbose = true );
}

//=============================================================================
StatusCode TupleToolApplyIsolationVetoDst::initialize() {
  if ( !TupleToolBase::initialize() ) return StatusCode::FAILURE;

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc() );
  if ( m_dva == nullptr )
    return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE );

  m_dist = tool<IDistanceCalculator>( "LoKi::DistanceCalculator", this );
  if ( !m_dist ) {
    Error( "Unable to retrieve the IDistanceCalculator tool" );
    return StatusCode::FAILURE;
  }
  m_pvReFitter = tool<IPVReFitter>( "AdaptivePVReFitter", this );
  m_pVertexFit = tool<IVertexFit>( "LoKi::VertexFitter", this );

  if ( !m_pVertexFit ) {
    Error( "Unable to retrieve the IVertexFit tool" );
    return StatusCode::FAILURE;
  }

  if ( m_trueID ) {
    std::vector<std::string> p2mcAssocTypes;
    p2mcAssocTypes.emplace_back( "DaVinciSmartAssociator" );
    p2mcAssocTypes.emplace_back( "MCMatchObjP2MCRelator" );
    for ( auto& p2mcAssocType : p2mcAssocTypes ) {
      m_p2mcAssocs.push_back(
          tool<IParticle2MCAssociator>( p2mcAssocType, this ) );
    }
  }

  m_Reader = new TMVA::Reader( "!Silent" );
  m_Reader->AddSpectator( "Track_TYPE", &type );
  m_Reader->AddVariable( "Track_MINIPCHI2", &minipchi2 );
  m_Reader->AddVariable( "Track_PT", &pt );
  m_Reader->AddVariable( "Track_OPENING", &opening );
  m_Reader->AddVariable( "Track_IPCHI2", &chi2 );
  m_Reader->AddVariable( "Track_FLIGHT", &newfdchi2 );
  m_Reader->AddVariable( "Track_DELTAFLIGHT", &deltafd );

  m_Reader->BookMVA( "BDT method", m_weightsName );
  if ( !m_Reader ) {
    Error( "Unable to retrieve the IVertexFit tool" );
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
StatusCode TupleToolApplyIsolationVetoDst::fill( const Particle*    mother,
                                                 const Particle*    P,
                                                 const std::string& head,
                                                 Tuples::Tuple&     tuple ) {
  const std::string prefix = fullName( head );
  Assert( P && mother && m_dist,
          "This should not happen, you are inside TupleToolVtxIsoln.cpp :(" );

  bool test = true;

  std::vector<const LHCb::Track*> daughtertracks;
  daughtertracks.clear();

  LHCb::Particle::ConstVector source;
  LHCb::Particle::ConstVector target;
  LHCb::Particle::ConstVector finalStates;
  LHCb::Particle::ConstVector parts2Vertex;
  LHCb::Particle::ConstVector parts2VertexDst;
  LHCb::Particle::ConstVector parts2VertexD;

  double                maxbdt[m_nWrite];
  const LHCb::Particle* SelParts[m_nWrite];

  vertexchi2 = P->endVertex()->chi2();
  parts2Vertex.clear();
  parts2VertexD.clear();
  parts2VertexDst.clear();

  const LHCb::Particle* theD = nullptr;

  LHCb::Particle::ConstVector bDaughters = mother->daughtersVector();
  // Mu loop
  for ( auto& bDaughter : bDaughters ) {
    if ( bDaughter->particleID().abspid() == 421 ) theD = bDaughter;
  }
  if ( !theD ) {
    warning() << "NOT FOUND D0 " << endmsg;
    tuple->column( prefix + "_ISOLATION_DstWindowOK", false );
    tuple->column( prefix + "_ISOLATION_DstWindowBDT", (float)0 );
    tuple->column( prefix + "_ISOLATION_DstWindowDELTAM", (float)0 );
    tuple->column( prefix + "_ISOLATION_DstWindowBDT2", (float)0 );
    tuple->column( prefix + "_ISOLATION_DstWindowDELTAM2", (float)0 );
    tuple->column( prefix + "_ISOLATION_DstWindowInWindow", (float)0 );

    return StatusCode( 0 );
  }
  parts2VertexDst.push_back( theD );
  for ( int i = 0; i < m_nWrite; i++ ) {
    maxbdt[i]   = -2;
    SelParts[i] = nullptr;
  }

  if ( P->isBasicParticle() ) {
    source.push_back( mother );
  } else {
    source.push_back( P );
  }
  LHCb::Vertex dv2;

  do {
    target.clear();
    for ( auto& isource : source ) {
      if ( !( isource->daughters().empty() ) ) {
        LHCb::Particle::ConstVector tmp = isource->daughtersVector();
        for ( auto& itmp : tmp ) {
          target.push_back( itmp );
          // Add the final states, i.e. particles with proto and ignoring
          // gammas
          if ( itmp->proto() && 22 != itmp->particleID().pid() ) {
            finalStates.push_back( itmp );
            daughtertracks.push_back( itmp->proto()->track() );
            if ( itmp->particleID().abspid() == 413 ) Dst_PT = itmp->pt();
          }
        }
      }  // if endVertex
    }    // isource
    source = target;
  } while ( !target.empty() );
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Final states size= " << finalStates.size() << endmsg;

  LHCb::Vertex v;

  if ( P->isBasicParticle() ) {
    parts2Vertex.push_back( P );
  } else {
    parts2Vertex  = P->daughtersVector();
    StatusCode sc = m_pVertexFit->fit( v, parts2Vertex );
  }

  // number below am IPCHI2 threshold isnt that useful, will probably remove it
  LHCb::Particle::ConstVector theParts;
  LHCb::Vertex                DstVertex;

  float DstBDT1 = -2;
  float DstBDT2 = -2;
  float DstM1   = -1;
  float DstM2   = -1;
  int   nWindow = 0;

  for ( auto& m_inputParticle : m_inputParticles ) {
    if ( !exist<LHCb::Particle::Range>( m_inputParticle + "/Particles" ) ) {
      if ( msgLevel( MSG::DEBUG ) )
        debug() << "No particles at " << m_inputParticle << " !!!!!" << endmsg;
      continue;
    }

    LHCb::Particle::Range parts =
        get<LHCb::Particle::Range>( m_inputParticle + "/Particles" );
    if ( msgLevel( MSG::DEBUG ) )
      debug() << "Getting particles from " << m_inputParticle << " with "
              << ( parts ).size() << " particles" << endmsg;
    for ( auto iparts : ( parts ) ) {
      const LHCb::Particle* part = iparts;

      if ( part->proto()->track()->type() < 5 &&
           !isTrackInDecay( part->proto()->track(), daughtertracks ) ) {
        ghostprob = part->proto()->track()->ghostProbability();
        if ( ghostprob > 0.5 ) {
          continue;
        }
        if ( part->proto()->track()->type() == 3 && opening <= 0.994 ) {
          continue;
        }
        if ( part->proto()->track()->type() == 4 && opening <= 0.98 ) {
          continue;
        }
        if ( part->proto()->track()->type() == 1 && opening <= 0.98 ) {
          continue;
        }

        LHCb::Vertex vtxWithExtraTrack;
        parts2Vertex.push_back( iparts );

        StatusCode sc3 = m_pVertexFit->fit( vtxWithExtraTrack, parts2Vertex );
        parts2Vertex.pop_back();

        opening   = getopening( part->proto()->track(), P );
        minipchi2 = getminipchi( part );
        newfdchi2 = getfdchi2( part->proto()->track(), vtxWithExtraTrack );
        oldfdchi2 = getfdchi2( part->proto()->track(), v );
        trackchi2 = part->proto()->track()->chi2PerDoF();
        deltafd   = log10( fabs( newfdchi2 - oldfdchi2 ) ) - 7;
        type      = part->proto()->track()->type();

        if ( newfdchi2 - oldfdchi2 < 0 ) deltafd = deltafd * -1.;
        newfdchi2 = log10( newfdchi2 );
        if ( part->proto()->track()->type() == 1 )
          pt = part->proto()->track()->momentum().z();
        else
          pt = part->proto()->track()->pt();

        float deltaMass = -1;
        if ( part->proto()->track()->type() == 3 ) {
          auto* newP = new LHCb::Particle( LHCb::ParticleID( -413 ) );
          parts2VertexDst.push_back( iparts );
          StatusCode sc4 =
              m_pVertexFit->fit( parts2VertexDst, DstVertex, *newP );
          parts2VertexDst.pop_back();
          deltaMass = newP->momentum().M() - theD->momentum().M();
        }

        bool deltaMassWindow = ( deltaMass > 140 && deltaMass < 160 );
        if ( deltaMassWindow ) {
          double     tmpip, tmpchi2;
          StatusCode dump =
              m_dist->distance( (const LHCb::Particle*)part,
                                (const LHCb::Vertex*)&v, tmpip, tmpchi2 );
          chi2 = tmpchi2;
          if ( chi2 > 100 ) continue;
          nWindow++;
          dummy        = 4000;
          float bdtval = m_Reader->EvaluateMVA( "BDT method" );
          if ( bdtval > DstBDT1 ) {
            DstBDT2 = DstBDT1;
            DstM2   = DstM1;
            DstBDT1 = bdtval;
            DstM1   = deltaMass;
          } else if ( bdtval > DstBDT2 ) {
            DstBDT2 = bdtval;
            DstM2   = deltaMass;
          }
          continue;
        }

        StatusCode sc = StatusCode::SUCCESS;
        double     tmpip, tmpchi2;
        StatusCode dump =
            m_dist->distance( (const LHCb::Particle*)part,
                              (const LHCb::Vertex*)&v, tmpip, tmpchi2 );
        chi2 = tmpchi2;

        if ( chi2 < 50 ) {
          dummy        = 4000;
          float bdtval = m_Reader->EvaluateMVA( "BDT method" );
          for ( int i = 0; i < m_nWrite; i++ ) {
            if ( bdtval > maxbdt[i] ) {
              for ( int j = m_nWrite; j > i; j-- ) {
                SelParts[j] = SelParts[j - 1];
                maxbdt[j]   = maxbdt[j - 1];
              }
              SelParts[i] = part;
              maxbdt[i]   = bdtval;
              break;
            }
          }
        }
      }
    }  // end particles loop
  }    // end particle types loop
  std::stringstream sout;
  for ( int i = 0; i < m_nWrite; i++ ) {
    std::string name;
    if ( i > 0 ) {
      sout << i + 1;
      name = sout.str();
    }
    writeParticle( SelParts[i], maxbdt[i], name, tuple, prefix, mother );
    sout.str( "" );
  }

  tuple->column( prefix + "_ISOLATION_DstWindowOK", true );
  tuple->column( prefix + "_ISOLATION_DstWindowBDT", DstBDT1 );
  tuple->column( prefix + "_ISOLATION_DstWindowDELTAM", DstM1 );
  tuple->column( prefix + "_ISOLATION_DstWindowBDT2", DstBDT2 );
  tuple->column( prefix + "_ISOLATION_DstWindowDELTAM2", DstM2 );
  tuple->column( prefix + "_ISOLATION_DstWindowInWindow", nWindow );

  return StatusCode( test );
}

//=========================================================================
const Vertex* TupleToolApplyIsolationVetoDst::originVertex(
    const Particle* top, const Particle* P ) const {
  if ( top == P || P->isBasicParticle() ) return nullptr;

  const SmartRefVector<LHCb::Particle>& dau = top->daughters();
  if ( dau.empty() ) return nullptr;

  SmartRefVector<LHCb::Particle>::const_iterator it;
  for ( it = dau.begin(); dau.end() != it; ++it ) {
    if ( P == *it ) return top->endVertex();  // I found the daughter
  }

  // vertex not yet found, get deeper in the decay:
  for ( it = dau.begin(); dau.end() != it; ++it ) {
    if ( P != *it && !( *it )->isBasicParticle() ) {
      const Vertex* vv = originVertex( *it, P );
      if ( vv ) return vv;
    }
  }

  return nullptr;
}

//=============================================================================
// Check if the track is already in the decay
//=============================================================================
bool TupleToolApplyIsolationVetoDst::isTrackInDecay(
    const LHCb::Track*                     track,
    const std::vector<const LHCb::Track*>& daughters ) {
  bool isInDecay = false;
  // loop over daughters
  for ( auto itrack : daughters ) {
    if ( itrack ) {
      if ( itrack == track ) {
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Track is in decay, skipping it" << endmsg;
        isInDecay = true;
      }
    }
  }  // end daughter loop

  return isInDecay;
}

//=============================================================================
// MINIPCHI2 for a track
//=============================================================================
double TupleToolApplyIsolationVetoDst::getminipchi(
    const LHCb::Particle* track ) {
  double                 minchi2 = -1;
  const RecVertex::Range PV      = m_dva->primaryVertices();

  if ( !PV.empty() ) {
    for ( auto pv : PV ) {
      double     ip, chi2;
      StatusCode test2 =
          m_dist->distance( (const LHCb::Particle*)track, pv, ip, chi2 );
      if ( ( chi2 < minchi2 ) || ( minchi2 < 0. ) ) {
        LHCb::RecVertex  newPV( *pv );
        StatusCode       scfit    = m_pvReFitter->remove( track, &newPV );
        LHCb::RecVertex* newPVPtr = (LHCb::RecVertex*)&newPV;
        test2                     = m_dist->distance( (LHCb::Particle*)track,
                                  (LHCb::VertexBase*)newPVPtr, ip, chi2 );
        minchi2                   = chi2;
      }
    }
  }

  return minchi2;
}

double TupleToolApplyIsolationVetoDst::getfdchi2( const LHCb::Track* track,
                                                  LHCb::Vertex       Vtx ) {
  double                 minchi2 = -1;
  double                 fdchi2  = -1;
  double                 fd;
  const RecVertex::Range PV = m_dva->primaryVertices();
  if ( !PV.empty() ) {
    for ( auto pv : PV ) {
      double     ip, chi2;
      StatusCode test2 =
          m_dist->distance( (const LHCb::Track*)track, pv, ip, chi2 );
      if ( ( chi2 < minchi2 ) || ( minchi2 < 0. ) ) {
        minchi2          = chi2;
        StatusCode test2 = m_dist->distance( pv, &Vtx, fd, fdchi2 );
      }
    }
  }

  return fdchi2;
}

//=============================================================================
// Opening angle for a track and particle
//=============================================================================
double TupleToolApplyIsolationVetoDst::getopening( const LHCb::Track*    track,
                                                   const LHCb::Particle* P ) {
  Gaudi::XYZVector A          = P->momentum().Vect();
  Gaudi::XYZVector B          = track->momentum();
  double           cosopening = A.Dot( B ) / std::sqrt( A.Mag2() * B.Mag2() );
  return cosopening;
}

//=============================================================================
// Store information in tuple for selected particles
//=============================================================================
void TupleToolApplyIsolationVetoDst::writeParticle(
    const LHCb::Particle* P, double bdt, const std::string& name,
    Tuples::Tuple& tuple, std::string prefix, const LHCb::Particle* Mother ) {
  if ( bdt > -1 ) {
    if ( m_trueID ) {
      const LHCb::MCParticle* mcp( nullptr );
      const LHCb::MCParticle* mother_mcp( nullptr );
      bool                    foundMCP        = false;
      bool                    foundMCP_mother = false;

      if ( P ) {
        // assignedPid = P->particleID().pid();
        if ( msgLevel( MSG::VERBOSE ) )
          verbose() << "Getting related MCP to " << P << endmsg;
        for ( auto& m_p2mcAssoc : m_p2mcAssocs ) {
          if ( !foundMCP ) mcp = m_p2mcAssoc->relatedMCP( P );
          if ( mcp ) {
            foundMCP = true;
          }
          if ( !foundMCP_mother )
            mother_mcp = m_p2mcAssoc->relatedMCP( Mother );
          if ( mother_mcp ) {
            foundMCP_mother = true;
          }
        }
        if ( msgLevel( MSG::VERBOSE ) )
          verbose() << "Got mcp " << mcp << endmsg;
      }
      float motherID = -1;
      if ( foundMCP )
        if ( mcp->mother() ) motherID = mcp->mother()->particleID().abspid();
      tuple->column( prefix + "_ISOLATION_MOTHERID" + name + m_outputSuffix,
                     motherID );
      float sameMother = 0;
      if ( foundMCP ) {
        if ( mother_mcp ) {
          const LHCb::MCParticle* head = mcp;
          for ( ;; ) {
            head = head->mother();
            if ( !head ) {
              break;
            }
            if ( ( head == mother_mcp ) ||
                 ( head->key() == mother_mcp->key() ) )
              sameMother = 1;
          }
        }
      }

      tuple->column( prefix + "_ISOLATION_SAMEMOTHER" + name + m_outputSuffix,
                     sameMother );
    }

  } else {
    if ( m_trueID ) {
      tuple->column( prefix + "_ISOLATION_MOTHERID" + name + m_outputSuffix,
                     (float)-1. );
      tuple->column( prefix + "_ISOLATION_SAMEMOTHER" + name + m_outputSuffix,
                     (float)-1 );
    }
  }
}
