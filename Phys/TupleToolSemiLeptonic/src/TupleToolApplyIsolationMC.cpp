// local
#include "TupleToolApplyIsolationMC.h"

// from Phys
#include "Kernel/GetIDVAlgorithm.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IDistanceCalculator.h"
#include "Kernel/IVertexFit.h"

// from Gaudi
#include "GaudiAlg/Tuple.h"

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

DECLARE_COMPONENT( TupleToolApplyIsolationMC )

using namespace LHCb;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolApplyIsolationMC::TupleToolApplyIsolationMC(
    const std::string& type, const std::string& name,
    const IInterface* parent )
    : TupleToolBase( type, name, parent ),
      m_dva( 0 ),
      m_dist( 0 ),
      m_pVertexFit( 0 ) {
  declareInterface<IParticleTupleTool>( this );

  m_inputParticles.emplace_back( "/Event/Phys/StdAllNoPIDsPions" );
  m_inputParticles.emplace_back( "/Event/Phys/StdNoPIDsUpPions" );
  m_inputParticles.emplace_back( "Phys/StdNoPIDsVeloPions" );

  m_p2mcAssocTypes.emplace_back( "DaVinciSmartAssociator" );
  m_p2mcAssocTypes.emplace_back( "MCMatchObjP2MCRelator" );

  // have not removed / added any of this yet
  declareProperty( "MaxDeltaChi2", m_deltaChi2 = 9.0 );
  declareProperty( "MaxChi2", m_Chi2 = 9.0 );
  declareProperty( "VertexFit", m_typeVertexFit = "default" );
  declareProperty( "InputParticles", m_inputParticles );
  declareProperty( "OutputSuffix", m_outputSuffix = "" );
  declareProperty( "WeightsFile", m_weightsName = "weights.xml" );
  declareProperty( "IP2MCPAssociatorTypes", m_p2mcAssocTypes );
}

//=============================================================================
StatusCode TupleToolApplyIsolationMC::initialize() {
  if ( !TupleToolBase::initialize() ) return StatusCode::FAILURE;

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc() );
  if ( m_dva == nullptr )
    return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE );

  m_dist      = tool<IDistanceCalculator>( "LoKi::DistanceCalculator", this );
  m_p2mcAssoc = tool<IParticle2MCAssociator>( "DaVinciSmartAssociator", this );
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

  m_Reader = new TMVA::Reader( "!Silent" );
  m_Reader->AddSpectator( "Track_TYPE", &type );
  m_Reader->AddVariable( "Track_MINIPCHI2", &minipchi2 );
  m_Reader->AddVariable( "Track_PT", &pt );
  m_Reader->AddVariable( "Track_OPENING", &opening );
  m_Reader->AddVariable( "Track_IPCHI2", &chi2 );
  m_Reader->AddVariable( "Track_FLIGHT", &newfdchi2 );
  m_Reader->AddVariable( "Track_DELTAFLIGHT", &deltafd );
  // dummy variables
  // m_Reader->AddVariable( "Bplus_PT",&dummy);
  // m_Reader->AddVariable( "Dst_PT",&Dst_PT);
  // m_Reader->AddVariable( "Bplus_ENDVERTEX_CHI2",&vertexchi2);
  // m_Reader->AddVariable( "Dst_ENDVERTEX_CHI2",&dummy);
  // m_Reader->AddVariable( "Dst_FDCHI2_OWNPV",&dummy);
  // m_Reader->AddVariable( "log(1-D_DIRA_OWNPV)",&dummy);
  // m_Reader->AddVariable( "log(1-Bplus_DIRA_OWNPV)",&dummy);

  m_Reader->BookMVA( "BDT method", m_weightsName );
  if ( !m_Reader ) {
    Error( "Unable to retrieve the IVertexFit tool" );
    return StatusCode::FAILURE;
  }

  m_p2mcAssocs.clear();
  for ( auto& m_p2mcAssocType : m_p2mcAssocTypes ) {
    m_p2mcAssocs.push_back(
        tool<IParticle2MCAssociator>( m_p2mcAssocType, this ) );
  }
  if ( m_p2mcAssocs.empty() ) {
    Error( "No MC associators configured" );
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
StatusCode TupleToolApplyIsolationMC::fill( const Particle*    mother,
                                            const Particle*    P,
                                            const std::string& head,
                                            Tuples::Tuple&     tuple ) {
  const std::string prefix = fullName( head );
  Assert( P && mother && m_dist,
          "This should not happen, you are inside TupleToolVtxIsoln.cpp :(" );

  bool  test        = true;
  float charge      = 0;
  float type        = 0;
  float px          = 0;
  float py          = 0;
  float pz          = 0;
  float pe          = 0;
  float ptransverse = 0;
  float eta         = 0;
  float pidk        = 0;
  float pidp        = 0;
  float nnk         = 0;
  float nnpi        = 0;
  float nnp         = 0;
  float nng         = 0;
  float ismuon      = 0;
  int   truepid     = 0;

  std::vector<const LHCb::Track*> daughtertracks;
  daughtertracks.clear();

  LHCb::Particle::ConstVector source;
  LHCb::Particle::ConstVector target;
  LHCb::Particle::ConstVector finalStates;
  LHCb::Particle::ConstVector parts2Vertex;
  LHCb::Particle::ConstVector parts2VertexD;

  double                      maxbdt = -2;
  double                      bdt2   = -2;
  double                      bdt3   = -2;
  const LHCb::Particle*       maxpart;
  const LHCb::Particle*       part2;
  const LHCb::Particle*       part3;

  vertexchi2 = P->endVertex()->chi2();
  parts2Vertex.clear();
  parts2VertexD.clear();

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

        for ( LHCb::Particle::ConstVector::const_iterator itmp = tmp.begin();
              itmp != tmp.end(); itmp++ ) {
          target.push_back( *itmp );
          // Add the final states, i.e. particles with proto and ignoring
          // gammas
          if ( ( *itmp )->proto() && 22 != ( *itmp )->particleID().pid() ) {
            finalStates.push_back( *itmp );
            daughtertracks.push_back( ( *itmp )->proto()->track() );
            if ( ( *itmp )->particleID().abspid() == 413 )
              Dst_PT = ( *itmp )->pt();
          }
        }
      }  // if endVertex
    }    // isource
    source = target;
  } while ( !target.empty() );
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Final states size= " << finalStates.size() << std::endl;

  LHCb::Vertex v;

  if ( P->isBasicParticle() ) {
    parts2Vertex.push_back( P );
  } else {
    parts2Vertex  = P->daughtersVector();
    StatusCode sc = m_pVertexFit->fit( v, parts2Vertex );
  }

  // number below am IPCHI2 threshold isnt that useful, will probably remove it
  LHCb::Particle::ConstVector theParts;

  for ( auto& m_inputParticle : m_inputParticles ) {
    if ( !exist<LHCb::Particle::Range>( m_inputParticle + "/Particles" ) ) {
      if ( msgLevel( MSG::DEBUG ) )
        debug() << "No particles at " << m_inputParticle << " !!!!!"
                << std::endl;
      continue;
    }

    LHCb::Particle::Range parts =
        get<LHCb::Particle::Range>( m_inputParticle + "/Particles" );
    if ( msgLevel( MSG::DEBUG ) )
      debug() << "Getting particles from " << m_inputParticle << " with "
              << ( parts ).size() << " particles" << std::endl;
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
        newfdchi2 =
            fabs( getfdchi2( part->proto()->track(), vtxWithExtraTrack ) );
        oldfdchi2 = fabs( getfdchi2( part->proto()->track(), v ) );
        trackchi2 = part->proto()->track()->chi2PerDoF();
        deltafd   = log10( fabs( newfdchi2 - oldfdchi2 ) ) - 7;
        type      = part->proto()->track()->type();

        if ( newfdchi2 - oldfdchi2 < 0 ) deltafd = deltafd * -1.;
        newfdchi2 = log10( newfdchi2 );
        if ( part->proto()->track()->type() == 1 )
          pt = part->proto()->track()->momentum().z();
        else
          pt = part->proto()->track()->pt();

        StatusCode sc = StatusCode::SUCCESS;
        double     tmpip, tmpchi2;
        StatusCode dump =
            m_dist->distance( (const LHCb::Particle*)part,
                              (const LHCb::Vertex*)&v, tmpip, tmpchi2 );

        chi2 = tmpchi2;
        if ( chi2 < 50 ) {
          dummy        = 4000;
          float bdtval = m_Reader->EvaluateMVA( "BDT method" );
          if ( bdtval > maxbdt ) {
            bdt3    = bdt2;
            bdt2    = maxbdt;
            maxbdt  = bdtval;
            part3   = part2;
            part2   = maxpart;
            maxpart = part;
          } else if ( bdtval > bdt2 ) {
            bdt3  = bdt2;
            bdt2  = bdtval;
            part3 = part2;
            part2 = part;
          } else if ( bdtval > bdt3 ) {
            bdt3  = bdtval;
            part3 = part;
          }
        }
      }
    }  // end particles loop
  }    // end particle types loop

  if ( maxbdt > -1 ) {
    pe          = maxpart->momentum().E();
    px          = maxpart->momentum().Px();
    py          = maxpart->momentum().Py();
    pz          = maxpart->momentum().Pz();
    ptransverse = maxpart->momentum().Pt();
    eta         = maxpart->momentum().Eta();
    pidk = maxpart->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp = maxpart->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp  = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk  = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNk, -1 );
    nnpi = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1 );
    nng  = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1 );
    if ( maxpart->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = maxpart->proto()->track()->charge();
    }
    type                   = maxpart->proto()->track()->type();
    const MuonPID* muonPID = maxpart->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;
    truepid = 0;

    const LHCb::MCParticle* mcmaxpart( nullptr );
    for ( auto& m_p2mcAssoc : m_p2mcAssocs ) {
      mcmaxpart = m_p2mcAssoc->relatedMCP( maxpart );
      if ( mcmaxpart ) {
        truepid = mcmaxpart->particleID().pid();
        break;
      }
    }
  }

  tuple->column( prefix + "_ISOLATION_BDT" + m_outputSuffix, maxbdt );
  tuple->column( prefix + "_ISOLATION_CHARGE" + m_outputSuffix, charge );
  tuple->column( prefix + "_ISOLATION_Type" + m_outputSuffix, type );
  tuple->column( prefix + "_ISOLATION_PE" + m_outputSuffix, pe );
  tuple->column( prefix + "_ISOLATION_PX" + m_outputSuffix, px );
  tuple->column( prefix + "_ISOLATION_PY" + m_outputSuffix, py );
  tuple->column( prefix + "_ISOLATION_PZ" + m_outputSuffix, pz );
  tuple->column( prefix + "_ISOLATION_PT" + m_outputSuffix, ptransverse );
  tuple->column( prefix + "_ISOLATION_ETA" + m_outputSuffix, eta );
  tuple->column( prefix + "_ISOLATION_PIDK" + m_outputSuffix, pidk );
  tuple->column( prefix + "_ISOLATION_PIDp" + m_outputSuffix, pidp );
  tuple->column( prefix + "_ISOLATION_NNk" + m_outputSuffix, nnk );
  tuple->column( prefix + "_ISOLATION_NNpi" + m_outputSuffix, nnpi );
  tuple->column( prefix + "_ISOLATION_NNp" + m_outputSuffix, nnp );
  tuple->column( prefix + "_ISOLATION_IsMuon" + m_outputSuffix, ismuon );
  tuple->column( prefix + "_ISOLATION_NNghost" + m_outputSuffix, nng );
  tuple->column( prefix + "_ISOLATION_TruePID" + m_outputSuffix, truepid );

  if ( bdt2 > -1 ) {
    pe          = part2->momentum().E();
    px          = part2->momentum().Px();
    py          = part2->momentum().Py();
    pz          = part2->momentum().Pz();
    ptransverse = part2->momentum().Pt();
    eta         = part2->momentum().Eta();
    pidk        = part2->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp        = part2->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp         = part2->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk         = part2->proto()->info( LHCb::ProtoParticle::ProbNNk, -1000 );
    nnpi        = part2->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1000 );
    nng = part2->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1000 );
    if ( part2->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = part2->proto()->track()->charge();
    }
    type                   = part2->proto()->track()->type();
    const MuonPID* muonPID = part2->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;
    truepid = 0;

    const LHCb::MCParticle* mcpart2( nullptr );
    for ( auto& m_p2mcAssoc : m_p2mcAssocs ) {
      mcpart2 = m_p2mcAssoc->relatedMCP( part2 );
      if ( mcpart2 ) {
        truepid = mcpart2->particleID().pid();
        break;
      }
    }
  }

  tuple->column( prefix + "_ISOLATION_BDT2" + m_outputSuffix, bdt2 );
  tuple->column( prefix + "_ISOLATION_CHARGE2" + m_outputSuffix, charge );
  tuple->column( prefix + "_ISOLATION_Type2" + m_outputSuffix, type );
  tuple->column( prefix + "_ISOLATION_PE2" + m_outputSuffix, pe );
  tuple->column( prefix + "_ISOLATION_PX2" + m_outputSuffix, px );
  tuple->column( prefix + "_ISOLATION_PY2" + m_outputSuffix, py );
  tuple->column( prefix + "_ISOLATION_PZ2" + m_outputSuffix, pz );
  tuple->column( prefix + "_ISOLATION_PT2" + m_outputSuffix, ptransverse );
  tuple->column( prefix + "_ISOLATION_ETA2" + m_outputSuffix, eta );
  tuple->column( prefix + "_ISOLATION_PIDK2" + m_outputSuffix, pidk );
  tuple->column( prefix + "_ISOLATION_PIDp2" + m_outputSuffix, pidp );
  tuple->column( prefix + "_ISOLATION_NNk2" + m_outputSuffix, nnk );
  tuple->column( prefix + "_ISOLATION_NNpi2" + m_outputSuffix, nnpi );
  tuple->column( prefix + "_ISOLATION_NNp2" + m_outputSuffix, nnp );
  tuple->column( prefix + "_ISOLATION_IsMuon2" + m_outputSuffix, ismuon );
  tuple->column( prefix + "_ISOLATION_NNghost2" + m_outputSuffix, nng );
  tuple->column( prefix + "_ISOLATION_TruePID2" + m_outputSuffix, truepid );

  if ( bdt3 > -1 ) {
    pe          = part3->momentum().E();
    px          = part3->momentum().Px();
    py          = part3->momentum().Py();
    pz          = part3->momentum().Pz();
    ptransverse = part3->momentum().Pt();
    eta         = part3->momentum().Eta();
    pidk        = part3->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp        = part3->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp         = part3->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk         = part3->proto()->info( LHCb::ProtoParticle::ProbNNk, -1000 );
    nnpi        = part3->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1000 );
    nng = part3->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1000 );
    if ( part3->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = part3->proto()->track()->charge();
    }
    type                   = part3->proto()->track()->type();
    const MuonPID* muonPID = part3->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;
    truepid = 0;

    const LHCb::MCParticle* mcpart3( nullptr );
    for ( auto& m_p2mcAssoc : m_p2mcAssocs ) {
      mcpart3 = m_p2mcAssoc->relatedMCP( part3 );
      if ( mcpart3 ) {
        truepid = mcpart3->particleID().pid();
        break;
      }
    }
  }

  tuple->column( prefix + "_ISOLATION_BDT3" + m_outputSuffix, bdt3 );
  tuple->column( prefix + "_ISOLATION_CHARGE3" + m_outputSuffix, charge );
  tuple->column( prefix + "_ISOLATION_Type3" + m_outputSuffix, type );
  tuple->column( prefix + "_ISOLATION_PE3" + m_outputSuffix, pe );
  tuple->column( prefix + "_ISOLATION_PX3" + m_outputSuffix, px );
  tuple->column( prefix + "_ISOLATION_PY3" + m_outputSuffix, py );
  tuple->column( prefix + "_ISOLATION_PZ3" + m_outputSuffix, pz );
  tuple->column( prefix + "_ISOLATION_PT3" + m_outputSuffix, ptransverse );
  tuple->column( prefix + "_ISOLATION_ETA3" + m_outputSuffix, eta );
  tuple->column( prefix + "_ISOLATION_PIDK3" + m_outputSuffix, pidk );
  tuple->column( prefix + "_ISOLATION_PIDp3" + m_outputSuffix, pidp );
  tuple->column( prefix + "_ISOLATION_NNk3" + m_outputSuffix, nnk );
  tuple->column( prefix + "_ISOLATION_NNpi3" + m_outputSuffix, nnpi );
  tuple->column( prefix + "_ISOLATION_NNp3" + m_outputSuffix, nnp );
  tuple->column( prefix + "_ISOLATION_IsMuon3" + m_outputSuffix, ismuon );
  tuple->column( prefix + "_ISOLATION_NNghost3" + m_outputSuffix, nng );
  tuple->column( prefix + "_ISOLATION_OldFDCHI2" + m_outputSuffix, oldfdchi2 );
  tuple->column( prefix + "_ISOLATION_TruePID3" + m_outputSuffix, truepid );

  return StatusCode( test );
}

//=========================================================================
const Vertex* TupleToolApplyIsolationMC::originVertex(
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
bool TupleToolApplyIsolationMC::isTrackInDecay(
    const LHCb::Track*                    track,
    const std::vector<const LHCb::Track*> daughters ) {
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
double TupleToolApplyIsolationMC::getminipchi( const LHCb::Particle* track ) {
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

double TupleToolApplyIsolationMC::getfdchi2( const LHCb::Track* track,
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
double TupleToolApplyIsolationMC::getopening( const LHCb::Track*    track,
                                              const LHCb::Particle* P ) {
  Gaudi::XYZVector A          = P->momentum().Vect();
  Gaudi::XYZVector B          = track->momentum();
  double           cosopening = A.Dot( B ) / std::sqrt( A.Mag2() * B.Mag2() );
  return cosopening;
}
