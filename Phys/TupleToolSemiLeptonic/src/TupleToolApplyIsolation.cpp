// local
#include "TupleToolApplyIsolation.h"

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

#include <cmath>
#include <functional>
#include <string>

DECLARE_COMPONENT( TupleToolApplyIsolation )

using namespace LHCb;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolApplyIsolation::TupleToolApplyIsolation( const std::string& type,
                                                  const std::string& name,
                                                  const IInterface*  parent )
    : TupleToolBase( type, name, parent ),
      m_dva( nullptr ),
      m_dist( nullptr ),
      m_pVertexFit( nullptr ) {
  declareInterface<IParticleTupleTool>( this );

  m_inputParticles.emplace_back( "/Event/Phys/StdAllNoPIDsPions" );
  m_inputParticles.emplace_back( "/Event/Phys/StdNoPIDsUpPions" );
  m_inputParticles.emplace_back( "Phys/StdNoPIDsVeloPions" );

  // have not removed / added any of this yet
  declareProperty( "MaxDeltaChi2", m_deltaChi2 = 9.0 );
  declareProperty( "MaxChi2", m_Chi2 = 9.0 );
  declareProperty( "VertexFit", m_typeVertexFit = "default" );
  declareProperty( "InputParticles", m_inputParticles );
  declareProperty( "OutputSuffix", m_outputSuffix = "" );
  declareProperty( "WeightsFile", m_weightsName = "weights.xml" );
}

//=============================================================================
StatusCode TupleToolApplyIsolation::initialize() {
  if ( !TupleToolBase::initialize() ) return StatusCode::FAILURE;

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc() );
  if ( m_dva == nullptr )
    return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE );

  m_dist      = tool<IDistanceCalculator>( "LoKi::DistanceCalculator", this );
  m_p2mcAssoc = tool<IParticle2MCAssociator>( "DaVinciSmartAssociator", this );
  if ( !m_dist ) {
    Error( "Unable to retrieve the IDistanceCalculator tool" ).ignore();
    return StatusCode::FAILURE;
  }

  m_pvReFitter = tool<IPVReFitter>( "AdaptivePVReFitter", this );
  m_pVertexFit = tool<IVertexFit>( "LoKi::VertexFitter", this );
  if ( !m_pVertexFit ) {
    Error( "Unable to retrieve the IVertexFit tool" ).ignore();
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
    Error( "Unable to retrieve the IVertexFit tool" ).ignore();
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
StatusCode TupleToolApplyIsolation::fill( const Particle*    mother,
                                          const Particle*    P,
                                          const std::string& head,
                                          Tuples::Tuple&     tuple ) {
  const std::string prefix = fullName( head );
  Assert( P && mother && m_dist,
          "This should not happen, you are inside TupleToolVtxIsoln.cpp :(" );

  bool  test   = true;
  float charge = 0;
  float type   = 0;
  float px     = 0;
  float py     = 0;
  float pz     = 0;
  float pe     = 0;
  float pidk   = 0;
  float pidp   = 0;
  float nnk    = 0;
  float nnpi   = 0;
  float nnp    = 0;
  float nng    = 0;
  float gprob  = 0;
  float ismuon = 0;

  std::vector<const LHCb::Track*> daughtertracks;
  daughtertracks.clear();

  LHCb::Particle::ConstVector source;
  LHCb::Particle::ConstVector target;
  LHCb::Particle::ConstVector finalStates;
  LHCb::Particle::ConstVector parts2Vertex;
  LHCb::Particle::ConstVector parts2VertexD;

  double                angle  = -2;
  double                angle2 = -2;
  double                angle3 = -2;
  double                angle4 = -2;
  double                angle5 = -2;
  double                maxchi2 = -99;
  double                mchi22  = -99;
  double                mchi23  = -99;
  double                mchi24  = -99;
  double                mchi25  = -99;
  double                maxbdt  = -2;
  double                bdt2    = -2;
  double                bdt3    = -2;
  double                bdt4    = -2;
  double                bdt5    = -2;
  int                   trueID  = 0;
  int                   _sc     = -999;
  int                   _sc2    = -999;
  int                   _sc3    = -999;
  int                   _sc4    = -999;
  int                   _sc5    = -999;
  const LHCb::Particle* maxpart = nullptr;
  const LHCb::Particle* part2 = nullptr;
  const LHCb::Particle* part3 = nullptr;
  const LHCb::Particle* part4 = nullptr;
  const LHCb::Particle* part5 = nullptr;

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
        opening = getopening( part->proto()->track(), P );
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
            bdt5    = bdt4;
            bdt4    = bdt3;
            bdt3    = bdt2;
            bdt2    = maxbdt;
            maxbdt  = bdtval;
            part5   = part4;
            part4   = part3;
            part3   = part2;
            part2   = maxpart;
            maxpart = part;
            mchi25  = mchi24;
            mchi24  = mchi23;
            mchi23  = mchi22;
            mchi22  = maxchi2;
            maxchi2 = tmpchi2;
            _sc5    = _sc4;
            _sc4    = _sc3;
            _sc3    = _sc2;
            _sc2    = _sc;
            _sc     = sc3.getCode();
            angle5  = angle4;
            angle4  = angle3;
            angle3  = angle2;
            angle2  = angle;
            angle   = opening;
          } else if ( bdtval > bdt2 ) {
            bdt5   = bdt4;
            bdt4   = bdt3;
            bdt3   = bdt2;
            bdt2   = bdtval;
            part5  = part4;
            part4  = part3;
            part3  = part2;
            part2  = part;
            mchi25 = mchi24;
            mchi24 = mchi23;
            mchi23 = mchi22;
            mchi22 = tmpchi2;
            _sc5   = _sc4;
            _sc4   = _sc3;
            _sc3   = _sc2;
            _sc2   = sc3.getCode();
            angle5 = angle4;
            angle4 = angle3;
            angle3 = angle2;
            angle2 = opening;
          } else if ( bdtval > bdt3 ) {
            bdt5   = bdt4;
            bdt4   = bdt3;
            bdt3   = bdtval;
            part5  = part4;
            part4  = part3;
            part3  = part;
            mchi25 = mchi24;
            mchi24 = mchi23;
            mchi23 = tmpchi2;
            _sc5   = _sc4;
            _sc4   = _sc3;
            _sc3   = sc3.getCode();
            angle5 = angle4;
            angle4 = angle3;
            angle3 = opening;
          } else if ( bdtval > bdt4 ) {
            bdt5   = bdt4;
            bdt4   = bdtval;
            part5  = part4;
            part4  = part;
            mchi25 = mchi24;
            mchi24 = tmpchi2;
            _sc5   = _sc4;
            _sc4   = sc3.getCode();
            angle5 = angle4;
            angle4 = opening;
          } else if ( bdtval > bdt5 ) {
            bdt5   = bdtval;
            part5  = part;
            mchi25 = tmpchi2;
            _sc5   = sc3.getCode();
            angle5 = opening;
          }
        }
      }
    }  // end particles loop
  }    // end particle types loop

  if ( maxbdt > -1 ) {
    pe    = maxpart->momentum().E();
    px    = maxpart->momentum().Px();
    py    = maxpart->momentum().Py();
    pz    = maxpart->momentum().Pz();
    pidk  = maxpart->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp  = maxpart->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp   = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk   = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNk, -1 );
    nnpi  = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1 );
    nng   = maxpart->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1 );
    gprob = maxpart->proto()->track()->ghostProbability();
    if ( maxpart->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = maxpart->proto()->track()->charge();
    }
    type                   = maxpart->proto()->track()->type();
    const MuonPID* muonPID = maxpart->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;

    const LHCb::MCParticle* mcp( nullptr );
    if ( msgLevel( MSG::VERBOSE ) )
      verbose() << "Getting related MCP to " << maxpart << endmsg;
    mcp = m_p2mcAssoc->relatedMCP( maxpart );
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Got mcp " << mcp << endmsg;
    trueID = ( mcp ? mcp->particleID().pid() : 0 );
  }

  test &= tuple->column( prefix + "_ISOLATION_CHI2" + m_outputSuffix, maxchi2 );
  test &= tuple->column( prefix + "_ISOLATION_ANGLE" + m_outputSuffix, angle );
  test &= tuple->column( prefix + "_ISOLATION_SC" + m_outputSuffix, _sc );
  test &= tuple->column( prefix + "_ISOLATION_BDT" + m_outputSuffix, maxbdt );
  test &= tuple->column( prefix + "_ISOLATION_CHARGE" + m_outputSuffix, charge );
  test &= tuple->column( prefix + "_ISOLATION_Type" + m_outputSuffix, type );
  test &= tuple->column( prefix + "_ISOLATION_PE" + m_outputSuffix, pe );
  test &= tuple->column( prefix + "_ISOLATION_PX" + m_outputSuffix, px );
  test &= tuple->column( prefix + "_ISOLATION_PY" + m_outputSuffix, py );
  test &= tuple->column( prefix + "_ISOLATION_PZ" + m_outputSuffix, pz );
  test &= tuple->column( prefix + "_ISOLATION_PIDK" + m_outputSuffix, pidk );
  test &= tuple->column( prefix + "_ISOLATION_PIDp" + m_outputSuffix, pidp );
  test &= tuple->column( prefix + "_ISOLATION_NNk" + m_outputSuffix, nnk );
  test &= tuple->column( prefix + "_ISOLATION_NNpi" + m_outputSuffix, nnpi );
  test &= tuple->column( prefix + "_ISOLATION_NNp" + m_outputSuffix, nnp );
  test &= tuple->column( prefix + "_ISOLATION_IsMuon" + m_outputSuffix, ismuon );
  test &= tuple->column( prefix + "_ISOLATION_NNghost" + m_outputSuffix, nng );
  test &= tuple->column( prefix + "_ISOLATION_GhostProb" + m_outputSuffix, gprob );
  test &= tuple->column( prefix + "_ISOLATION_TRUEID" + m_outputSuffix, trueID );

  if ( bdt2 > -1 ) {
    pe    = part2->momentum().E();
    px    = part2->momentum().Px();
    py    = part2->momentum().Py();
    pz    = part2->momentum().Pz();
    pidk  = part2->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp  = part2->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp   = part2->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk   = part2->proto()->info( LHCb::ProtoParticle::ProbNNk, -1000 );
    nnpi  = part2->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1000 );
    nng   = part2->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1000 );
    gprob = part2->proto()->track()->ghostProbability();
    if ( part2->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = part2->proto()->track()->charge();
    }
    type                   = part2->proto()->track()->type();
    const MuonPID* muonPID = part2->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;

    const LHCb::MCParticle* mcp( nullptr );
    if ( msgLevel( MSG::VERBOSE ) )
      verbose() << "Getting related MCP to " << part2 << endmsg;
    mcp = m_p2mcAssoc->relatedMCP( part2 );
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Got mcp " << mcp << endmsg;
    trueID = ( mcp ? mcp->particleID().pid() : 0 );
  }

  test &= tuple->column( prefix + "_ISOLATION_CHI22" + m_outputSuffix, mchi22 );
  test &= tuple->column( prefix + "_ISOLATION_SC2" + m_outputSuffix, _sc2 );
  test &= tuple->column( prefix + "_ISOLATION_ANGLE2" + m_outputSuffix, angle2 );
  test &= tuple->column( prefix + "_ISOLATION_BDT2" + m_outputSuffix, bdt2 );
  test &= tuple->column( prefix + "_ISOLATION_CHARGE2" + m_outputSuffix, charge );
  test &= tuple->column( prefix + "_ISOLATION_Type2" + m_outputSuffix, type );
  test &= tuple->column( prefix + "_ISOLATION_PE2" + m_outputSuffix, pe );
  test &= tuple->column( prefix + "_ISOLATION_PX2" + m_outputSuffix, px );
  test &= tuple->column( prefix + "_ISOLATION_PY2" + m_outputSuffix, py );
  test &= tuple->column( prefix + "_ISOLATION_PZ2" + m_outputSuffix, pz );
  test &= tuple->column( prefix + "_ISOLATION_PIDK2" + m_outputSuffix, pidk );
  test &= tuple->column( prefix + "_ISOLATION_PIDp2" + m_outputSuffix, pidp );
  test &= tuple->column( prefix + "_ISOLATION_NNk2" + m_outputSuffix, nnk );
  test &= tuple->column( prefix + "_ISOLATION_NNpi2" + m_outputSuffix, nnpi );
  test &= tuple->column( prefix + "_ISOLATION_NNp2" + m_outputSuffix, nnp );
  test &= tuple->column( prefix + "_ISOLATION_IsMuon2" + m_outputSuffix, ismuon );
  test &= tuple->column( prefix + "_ISOLATION_NNghost2" + m_outputSuffix, nng );
  test &= tuple->column( prefix + "_ISOLATION_GhostProb2" + m_outputSuffix, gprob );
  test &= tuple->column( prefix + "_ISOLATION_TRUEID2" + m_outputSuffix, trueID );

  if ( bdt3 > -1 ) {
    pe    = part3->momentum().E();
    px    = part3->momentum().Px();
    py    = part3->momentum().Py();
    pz    = part3->momentum().Pz();
    pidk  = part3->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp  = part3->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp   = part3->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk   = part3->proto()->info( LHCb::ProtoParticle::ProbNNk, -1000 );
    nnpi  = part3->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1000 );
    nng   = part3->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1000 );
    gprob = part3->proto()->track()->ghostProbability();
    if ( part3->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = part3->proto()->track()->charge();
    }
    type                   = part3->proto()->track()->type();
    const MuonPID* muonPID = part3->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;

    const LHCb::MCParticle* mcp( nullptr );
    if ( msgLevel( MSG::VERBOSE ) )
      verbose() << "Getting related MCP to " << part3 << endmsg;
    mcp = m_p2mcAssoc->relatedMCP( part3 );
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Got mcp " << mcp << endmsg;
    trueID = ( mcp ? mcp->particleID().pid() : 0 );
  }

  test &= tuple->column( prefix + "_ISOLATION_CHI23" + m_outputSuffix, mchi23 );
  test &= tuple->column( prefix + "_ISOLATION_SC3" + m_outputSuffix, _sc3 );
  test &= tuple->column( prefix + "_ISOLATION_BDT3" + m_outputSuffix, bdt3 );
  test &= tuple->column( prefix + "_ISOLATION_ANGLE3" + m_outputSuffix, angle3 );
  test &= tuple->column( prefix + "_ISOLATION_CHARGE3" + m_outputSuffix, charge );
  test &= tuple->column( prefix + "_ISOLATION_Type3" + m_outputSuffix, type );
  test &= tuple->column( prefix + "_ISOLATION_PE3" + m_outputSuffix, pe );
  test &= tuple->column( prefix + "_ISOLATION_PX3" + m_outputSuffix, px );
  test &= tuple->column( prefix + "_ISOLATION_PY3" + m_outputSuffix, py );
  test &= tuple->column( prefix + "_ISOLATION_PZ3" + m_outputSuffix, pz );
  test &= tuple->column( prefix + "_ISOLATION_PIDK3" + m_outputSuffix, pidk );
  test &= tuple->column( prefix + "_ISOLATION_PIDp3" + m_outputSuffix, pidp );
  test &= tuple->column( prefix + "_ISOLATION_NNk3" + m_outputSuffix, nnk );
  test &= tuple->column( prefix + "_ISOLATION_NNpi3" + m_outputSuffix, nnpi );
  test &= tuple->column( prefix + "_ISOLATION_NNp3" + m_outputSuffix, nnp );
  test &= tuple->column( prefix + "_ISOLATION_IsMuon3" + m_outputSuffix, ismuon );
  test &= tuple->column( prefix + "_ISOLATION_NNghost3" + m_outputSuffix, nng );
  test &= tuple->column( prefix + "_ISOLATION_GhostProb3" + m_outputSuffix, gprob );
  test &= tuple->column( prefix + "_ISOLATION_TRUEID3" + m_outputSuffix, trueID );

  if ( bdt4 > -1 ) {
    pe   = part4->momentum().E();
    px   = part4->momentum().Px();
    py   = part4->momentum().Py();
    pz   = part4->momentum().Pz();
    pidk = part4->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp = part4->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp  = part4->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk  = part4->proto()->info( LHCb::ProtoParticle::ProbNNk, -1000 );
    nnpi = part4->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1000 );
    nng  = part4->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1000 );
    gprob = part4->proto()->track()->ghostProbability();
    if ( part4->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = part4->proto()->track()->charge();
    }
    type                   = part4->proto()->track()->type();
    const MuonPID* muonPID = part4->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;

    const LHCb::MCParticle* mcp( nullptr );
    if ( msgLevel( MSG::VERBOSE ) )
      verbose() << "Getting related MCP to " << part4 << endmsg;
    mcp = m_p2mcAssoc->relatedMCP( part4 );
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Got mcp " << mcp << endmsg;
    trueID = ( mcp ? mcp->particleID().pid() : 0 );
  }

  test &= tuple->column( prefix + "_ISOLATION_CHI24" + m_outputSuffix, mchi24 );
  test &= tuple->column( prefix + "_ISOLATION_SC4" + m_outputSuffix, _sc4 );
  test &= tuple->column( prefix + "_ISOLATION_BDT4" + m_outputSuffix, bdt4 );
  test &= tuple->column( prefix + "_ISOLATION_ANGLE4" + m_outputSuffix, angle4 );
  test &= tuple->column( prefix + "_ISOLATION_CHARGE4" + m_outputSuffix, charge );
  test &= tuple->column( prefix + "_ISOLATION_Type4" + m_outputSuffix, type );
  test &= tuple->column( prefix + "_ISOLATION_PE4" + m_outputSuffix, pe );
  test &= tuple->column( prefix + "_ISOLATION_PX4" + m_outputSuffix, px );
  test &= tuple->column( prefix + "_ISOLATION_PY4" + m_outputSuffix, py );
  test &= tuple->column( prefix + "_ISOLATION_PZ4" + m_outputSuffix, pz );
  test &= tuple->column( prefix + "_ISOLATION_PIDK4" + m_outputSuffix, pidk );
  test &= tuple->column( prefix + "_ISOLATION_PIDp4" + m_outputSuffix, pidp );
  test &= tuple->column( prefix + "_ISOLATION_NNk4" + m_outputSuffix, nnk );
  test &= tuple->column( prefix + "_ISOLATION_NNpi4" + m_outputSuffix, nnpi );
  test &= tuple->column( prefix + "_ISOLATION_NNp4" + m_outputSuffix, nnp );
  test &= tuple->column( prefix + "_ISOLATION_IsMuon4" + m_outputSuffix, ismuon );
  test &= tuple->column( prefix + "_ISOLATION_NNghost4" + m_outputSuffix, nng );
  test &= tuple->column( prefix + "_ISOLATION_GhostProb4" + m_outputSuffix, gprob );
  test &= tuple->column( prefix + "_ISOLATION_TRUEID4" + m_outputSuffix, trueID );

  if ( bdt5 > -1 ) {
    pe    = part5->momentum().E();
    px    = part5->momentum().Px();
    py    = part5->momentum().Py();
    pz    = part5->momentum().Pz();
    pidk  = part5->proto()->info( LHCb::ProtoParticle::CombDLLk, -1000 );
    pidp  = part5->proto()->info( LHCb::ProtoParticle::CombDLLp, -1000 );
    nnp   = part5->proto()->info( LHCb::ProtoParticle::ProbNNp, -1 );
    nnk   = part5->proto()->info( LHCb::ProtoParticle::ProbNNk, -1000 );
    nnpi  = part5->proto()->info( LHCb::ProtoParticle::ProbNNpi, -1000 );
    nng   = part5->proto()->info( LHCb::ProtoParticle::ProbNNghost, -1000 );
    gprob = part5->proto()->track()->ghostProbability();
    if ( part5->proto()->track()->type() == 1 ) {
      charge = 0;
    } else {
      charge = part5->proto()->track()->charge();
    }
    type                   = part5->proto()->track()->type();
    const MuonPID* muonPID = part5->proto()->muonPID();
    ismuon                 = muonPID ? muonPID->IsMuon() : false;

    const LHCb::MCParticle* mcp( nullptr );
    if ( msgLevel( MSG::VERBOSE ) )
      verbose() << "Getting related MCP to " << part5 << endmsg;
    mcp = m_p2mcAssoc->relatedMCP( part5 );
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Got mcp " << mcp << endmsg;
    trueID = ( mcp ? mcp->particleID().pid() : 0 );
  }

  test &= tuple->column( prefix + "_ISOLATION_CHI25" + m_outputSuffix, mchi25 );
  test &= tuple->column( prefix + "_ISOLATION_SC5" + m_outputSuffix, _sc5 );
  test &= tuple->column( prefix + "_ISOLATION_BDT5" + m_outputSuffix, bdt5 );
  test &= tuple->column( prefix + "_ISOLATION_ANGLE5" + m_outputSuffix, angle5 );
  test &= tuple->column( prefix + "_ISOLATION_CHARGE5" + m_outputSuffix, charge );
  test &= tuple->column( prefix + "_ISOLATION_Type5" + m_outputSuffix, type );
  test &= tuple->column( prefix + "_ISOLATION_PE5" + m_outputSuffix, pe );
  test &= tuple->column( prefix + "_ISOLATION_PX5" + m_outputSuffix, px );
  test &= tuple->column( prefix + "_ISOLATION_PY5" + m_outputSuffix, py );
  test &= tuple->column( prefix + "_ISOLATION_PZ5" + m_outputSuffix, pz );
  test &= tuple->column( prefix + "_ISOLATION_PIDK5" + m_outputSuffix, pidk );
  test &= tuple->column( prefix + "_ISOLATION_PIDp5" + m_outputSuffix, pidp );
  test &= tuple->column( prefix + "_ISOLATION_NNk5" + m_outputSuffix, nnk );
  test &= tuple->column( prefix + "_ISOLATION_NNpi5" + m_outputSuffix, nnpi );
  test &= tuple->column( prefix + "_ISOLATION_NNp5" + m_outputSuffix, nnp );
  test &= tuple->column( prefix + "_ISOLATION_IsMuon5" + m_outputSuffix, ismuon );
  test &= tuple->column( prefix + "_ISOLATION_NNghost5" + m_outputSuffix, nng );
  test &= tuple->column( prefix + "_ISOLATION_GhostProb5" + m_outputSuffix, gprob );
  test &= tuple->column( prefix + "_ISOLATION_TRUEID5" + m_outputSuffix, trueID );

  return StatusCode( test );
}

//=========================================================================
const Vertex* TupleToolApplyIsolation::originVertex(
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
bool TupleToolApplyIsolation::isTrackInDecay(
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
  }

  return isInDecay;
}

//=============================================================================
// MINIPCHI2 for a track
//=============================================================================
double TupleToolApplyIsolation::getminipchi( const LHCb::Particle* track ) {
  double                 minchi2 = -1;
  const RecVertex::Range PV      = m_dva->primaryVertices();

  if ( !PV.empty() ) {
    for ( auto pv : PV ) {
      double     ip, chi2;
      StatusCode test2 =
          m_dist->distance( (const LHCb::Particle*)track, pv, ip, chi2 );
      if ( ( chi2 < minchi2 ) || ( minchi2 < 0. ) ) {
        LHCb::RecVertex newPV( *pv );
        StatusCode      scfit    = m_pvReFitter->remove( track, &newPV );
        auto*           newPVPtr = (LHCb::RecVertex*)&newPV;
        test2                    = m_dist->distance( (LHCb::Particle*)track,
                                  (LHCb::VertexBase*)newPVPtr, ip, chi2 );
        minchi2                  = chi2;
      }
    }
  }

  return minchi2;
}

double TupleToolApplyIsolation::getfdchi2( const LHCb::Track* track,
                                           LHCb::Vertex       Vtx ) {
  double                 minchi2 = -1;
  double                 fdchi2  = -1;
  double                 fd;
  const RecVertex::Range PV = m_dva->primaryVertices();

  if ( !PV.empty() ) {
    for ( auto pv : PV ) {
      double     ip = 0, chi2;
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
double TupleToolApplyIsolation::getopening( const LHCb::Track*    track,
                                            const LHCb::Particle* P ) {
  Gaudi::XYZVector A          = P->momentum().Vect();
  Gaudi::XYZVector B          = track->momentum();
  double           cosopening = A.Dot( B ) / std::sqrt( A.Mag2() * B.Mag2() );
  return cosopening;
}
