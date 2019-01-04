// Include files

// from Gaudi
#include "GaudiKernel/ToolFactory.h"

// local
#include "TupleToolApplyIsolationVetoDst.h"

#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>
#include <Kernel/IVertexFit.h>

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/MCParticle.h"
#include "Event/Particle.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "Kernel/IPVReFitter.h"
#include "Linker/LinkerTable.h"
#include "TrackInterfaces/ITrackVertexer.h"
#include <functional>

#include <map>
#include <string>

#include "TFile.h"
#include "TPluginManager.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

//#include "TMVA/Reader.h"
//#include "TMVA/Config.h"
//#include "TMVA/Tools.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolVtxIsoln
//
// @author Mitesh Patel, Patrick Koppenburg
// @date   2008-04-15
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_TOOL_FACTORY(TupleToolApplyIsolationVetoDst);

using namespace LHCb;
//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolApplyIsolationVetoDst::TupleToolApplyIsolationVetoDst(const std::string& type,
    const std::string& name,
    const IInterface* parent)
    : TupleToolBase(type, name, parent)
    , m_dva(0)
    , m_dist(0)
    , m_pVertexFit(0)
{
    declareInterface<IParticleTupleTool>(this);

    m_inputParticles.push_back("/Event/Phys/StdAllNoPIDsPions");
    m_inputParticles.push_back("/Event/Phys/StdNoPIDsUpPions");
    m_inputParticles.push_back("Phys/StdNoPIDsVeloPions");
    //m_inputParticles.push_back("/Event/Phys/StdNoPIDsVeloElectrons");

    declareProperty("VertexFit", m_typeVertexFit = "default");
    declareProperty("InputParticles", m_inputParticles);
    declareProperty("OutputSuffix", m_outputSuffix = "");
    declareProperty("WeightsFile", m_weightsName = "weights.xml");
    declareProperty("NWrite", m_nWrite = 3);
    declareProperty("TrueIDs", m_trueID = false);
    declareProperty("Verbose", m_verbose = true);
}

//=============================================================================

StatusCode TupleToolApplyIsolationVetoDst::initialize()
{
    if (!TupleToolBase::initialize())
        return StatusCode::FAILURE;

    m_dva = Gaudi::Utils::getIDVAlgorithm(contextSvc());
    if (0 == m_dva)
        return Error("Couldn't get parent DVAlgorithm",
            StatusCode::FAILURE);
    m_dist = tool<IDistanceCalculator>("LoKi::DistanceCalculator", this);
    if (!m_dist) {
        Error("Unable to retrieve the IDistanceCalculator tool");
        return StatusCode::FAILURE;
    }
    m_pvReFitter = tool<IPVReFitter>("AdaptivePVReFitter", this);
    //m_pVertexFit= m_dva->vertexFitter();
    m_pVertexFit = tool<IVertexFit>("LoKi::VertexFitter", this);
    //m_pVertexFit= tool<ITrackVertexer>

    if (!m_pVertexFit) {
        Error("Unable to retrieve the IVertexFit tool");
        return StatusCode::FAILURE;
    }

    if (m_trueID) {

        std::vector<std::string> p2mcAssocTypes;
        p2mcAssocTypes.push_back("DaVinciSmartAssociator");
        p2mcAssocTypes.push_back("MCMatchObjP2MCRelator");
        for (std::vector<std::string>::const_iterator iMCAss = p2mcAssocTypes.begin();
             iMCAss != p2mcAssocTypes.end(); ++iMCAss) {
            m_p2mcAssocs.push_back(tool<IParticle2MCAssociator>(*iMCAss, this));
        }
    }

    m_Reader = new TMVA::Reader("!Silent");
    m_Reader->AddSpectator("Track_TYPE", &type);
    m_Reader->AddVariable("Track_MINIPCHI2", &minipchi2);
    m_Reader->AddVariable("Track_PT", &pt);
    m_Reader->AddVariable("Track_OPENING", &opening);
    m_Reader->AddVariable("Track_IPCHI2", &chi2);
    m_Reader->AddVariable("Track_FLIGHT", &newfdchi2);
    m_Reader->AddVariable("Track_DELTAFLIGHT", &deltafd);
    //dummy variables
    //m_Reader->AddVariable( "Bplus_PT",&dummy);
    //m_Reader->AddVariable( "Dst_PT",&Dst_PT);
    //m_Reader->AddVariable( "Bplus_ENDVERTEX_CHI2",&vertexchi2);
    //m_Reader->AddVariable( "Dst_ENDVERTEX_CHI2",&dummy);
    //m_Reader->AddVariable( "Dst_FDCHI2_OWNPV",&dummy);
    //m_Reader->AddVariable( "log(1-D_DIRA_OWNPV)",&dummy);
    //m_Reader->AddVariable( "log(1-Bplus_DIRA_OWNPV)",&dummy);

    //reader->AddVariable("Bplus_PT",&Bplus_PTf);
    //reader->AddVariable("Dst_PT",&Dst_PTf);
    m_Reader->BookMVA("BDT method", m_weightsName);

    if (!m_Reader) {
        Error("Unable to retrieve the IVertexFit tool");
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}

//=============================================================================

StatusCode TupleToolApplyIsolationVetoDst::fill(const Particle* mother, const Particle* P, const std::string& head, Tuples::Tuple& tuple)
{

    const std::string prefix = fullName(head);
    Assert(P && mother && m_dist, "This should not happen, you are inside TupleToolVtxIsoln.cpp :(");

    bool test = true;

    /*
  const LHCb::Vertex* vtx;
  if (P->isBasicParticle()){
    vtx = mother->endVertex(); 
  }
  else{
    vtx = P->endVertex();
    
  }
  debug()<<"vertex for P, ID " <<P->particleID().pid()<<" = " <<vtx<<" at "<<vtx->position()<<  endmsg;
  if( !vtx ){
    Error("Can't retrieve the  vertex for " + prefix );
    return StatusCode::FAILURE;
  }
  */
    std::vector<const LHCb::Track*> daughtertracks;
    daughtertracks.clear();

    LHCb::Particle::ConstVector source;
    LHCb::Particle::ConstVector target;
    LHCb::Particle::ConstVector finalStates;
    LHCb::Particle::ConstVector parts2Vertex;
    LHCb::Particle::ConstVector parts2VertexDst;
    LHCb::Particle::ConstVector parts2VertexD;
    double maxbdt[m_nWrite];
    const LHCb::Particle* SelParts[m_nWrite];
    vertexchi2 = P->endVertex()->chi2();
    parts2Vertex.clear();
    parts2VertexD.clear();
    parts2VertexDst.clear();

    const LHCb::Particle* theD = NULL;

    LHCb::Particle::ConstVector bDaughters = mother->daughtersVector();
    //Mu loop
    int iloop = 0;
    //      warning() << "Head PID " << mother->particleID().pid() << endmsg;
    for (LHCb::Particle::ConstVector::const_iterator itmp = bDaughters.begin();
         itmp != bDaughters.end(); ++itmp) {
        //       warning() << "PID " << (*itmp)->particleID().pid() << endmsg;
        //       warning() << "Mass " << (*itmp)->momentum().M() << endmsg;
        if ((*itmp)->particleID().abspid() == 421)
            theD = (*itmp);
    }
    if (!theD)
        warning() << "NOT FOUND D0 " << endmsg;
    parts2VertexDst.push_back(theD);
    for (int i = 0; i < m_nWrite; i++) {
        maxbdt[i] = -2;
        SelParts[i] = NULL;
    }

    //   const LHCb::Particle* prefix = P;
    if (P->isBasicParticle()) {
        source.push_back(mother);
    } else {
        source.push_back(P);
    }
    LHCb::Vertex dv2;

    do {
        target.clear();
        for (LHCb::Particle::ConstVector::const_iterator isource = source.begin();
             isource != source.end(); isource++) {

            if (!((*isource)->daughters().empty())) {

                LHCb::Particle::ConstVector tmp = (*isource)->daughtersVector();

                for (LHCb::Particle::ConstVector::const_iterator itmp = tmp.begin();
                     itmp != tmp.end(); itmp++) {
                    target.push_back(*itmp);
                    // Add the final states, i.e. particles with proto and ignoring gammas
                    if ((*itmp)->proto() && 22 != (*itmp)->particleID().pid()) {
                        finalStates.push_back(*itmp);
                        daughtertracks.push_back((*itmp)->proto()->track());
                        if ((*itmp)->particleID().abspid() == 413)
                            Dst_PT = (*itmp)->pt();
                    }
                }
            } // if endVertex
        } // isource
        source = target;
    } while (target.size() > 0);
    if (msgLevel(MSG::DEBUG))
        debug() << "Final states size= " << finalStates.size() << endmsg;
    //warning() << " D VERTEX CHI2 " << dv2.chi2() << " NDOF " << dv2.nDoF() << endmsg;
    //warning() << "DAUGHTER SIZE " << daughtertracks.size() << endmsg;
    //default gives best tracks (why would you default to anything less than the best?)

    LHCb::Vertex v;
    //double chi2ndof = 0;//oldvtx->chi2();
    //int ndof = 0;//oldvtx->nDoF();

    if (P->isBasicParticle()) {
        parts2Vertex.push_back(P);
    } else {
        parts2Vertex = P->daughtersVector();
        StatusCode sc = m_pVertexFit->fit(v, parts2Vertex);
    }

    //hack due to lack of programming skill
    //v=new LHCb::Vertex(v2);
    //********Loop through tracks************

    //number below am IPCHI2 threshold isnt that useful, will probably remove it

    LHCb::Particle::ConstVector theParts;

    LHCb::Vertex DstVertex;

    float DstBDT1 = -2;
    float DstBDT2 = -2;
    float DstM1 = -1;
    float DstM2 = -1;
    int nWindow = 0;

    for (std::vector<std::string>::iterator i = m_inputParticles.begin();
         i != m_inputParticles.end(); ++i) {

        if (!exist<LHCb::Particle::Range>(*i + "/Particles")) {
            if (msgLevel(MSG::DEBUG))
                debug() << "No particles at " << *i << " !!!!!" << endmsg;
            continue;
        }

        LHCb::Particle::Range parts = get<LHCb::Particle::Range>(*i + "/Particles");
        if (msgLevel(MSG::DEBUG))
            debug() << "Getting particles from " << *i
                    << " with " << (parts).size() << " particles" << endmsg;
        //warning() << "Getting particles from " << *i
        //                                  << " with " << (parts).size() << " particles" << endmsg;
        for (LHCb::Particle::Range::const_iterator iparts = (parts).begin();
             iparts != (parts).end(); ++iparts) {
            const LHCb::Particle* part = (*iparts);

            //if(isTrackInDecay(part->proto()->track(),daughtertracks)) warning() << "FOUND DAUGHTER TRACK" << endmsg;
            if (part->proto()->track()->type() < 5 && !isTrackInDecay(part->proto()->track(), daughtertracks)) {
                LHCb::Vertex vtxWithExtraTrack;
                parts2Vertex.push_back(*iparts);
                StatusCode sc3 = m_pVertexFit->fit(vtxWithExtraTrack, parts2Vertex);
                parts2Vertex.pop_back();
                opening = getopening(part->proto()->track(), P);
                minipchi2 = getminipchi(part);
                newfdchi2 = getfdchi2(part->proto()->track(), vtxWithExtraTrack);
                oldfdchi2 = getfdchi2(part->proto()->track(), v);
                ghostprob = part->proto()->track()->ghostProbability();
                trackchi2 = part->proto()->track()->chi2PerDoF();
                deltafd = log10(fabs(newfdchi2 - oldfdchi2)) - 7;
                type = part->proto()->track()->type();
                if (newfdchi2 - oldfdchi2 < 0)
                    deltafd = deltafd * -1.;
                //warning() << "DELTAFD " << deltafd << endmsg;
                newfdchi2 = log10(newfdchi2);
                if (part->proto()->track()->type() == 1)
                    pt = part->proto()->track()->momentum().z();
                else
                    pt = part->proto()->track()->pt();

                float deltaMass = -1;
                if (part->proto()->track()->type() == 3) {
                    LHCb::Particle* newP = new LHCb::Particle(LHCb::ParticleID(-413));
                    parts2VertexDst.push_back(*iparts);
                    StatusCode sc4 = m_pVertexFit->fit(parts2VertexDst, DstVertex, *newP);
                    parts2VertexDst.pop_back();
                    deltaMass = newP->momentum().M() - theD->momentum().M();
                    //       warning() << "Dstar mass"<< newP->momentum().M() << " D mass " << theD->momentum().M() <<  "deltaMass " << deltaMass << " fit stats " <<  sc4 << endmsg;
                }
                //     warning() << "deltaMass " << deltaMass << " track type " << part->proto()->track()->type() << endr
                bool deltaMassWindow = (deltaMass > 140 && deltaMass < 160);

                if (deltaMassWindow) {
                    double tmpip, tmpchi2;
                    StatusCode dump = m_dist->distance((const LHCb::Particle*)part, (const LHCb::Vertex*)&v, tmpip, tmpchi2);
                    chi2 = tmpchi2;
                    if (chi2 > 100)
                        continue;
                    nWindow++;
                    dummy = 4000;
                    float bdtval = m_Reader->EvaluateMVA("BDT method");
                    if (bdtval > DstBDT1) {
                        DstBDT2 = DstBDT1;
                        DstM2 = DstM1;
                        DstBDT1 = bdtval;
                        DstM1 = deltaMass;
                    } else if (bdtval > DstBDT2) {
                        DstBDT2 = bdtval;
                        DstM2 = deltaMass;
                    }
                    continue;
                }

                //     if(!deltaMassWindow){
                if (ghostprob > 0.5) {
                    continue;
                }
                if (part->proto()->track()->type() == 3 && !(opening > 0.994)) {
                    continue;
                }
                if (part->proto()->track()->type() == 4 && !(opening > 0.98)) {
                    continue;
                }
                if (part->proto()->track()->type() == 1 && !(opening > 0.98)) {
                    continue;
                }
                //     }
                //if(track->info(LHCb::Track::CloneDist, -1.) > 0){continue;}
                StatusCode sc = StatusCode::SUCCESS;
                double tmpip, tmpchi2;
                StatusCode dump = m_dist->distance((const LHCb::Particle*)part, (const LHCb::Vertex*)&v, tmpip, tmpchi2);
                chi2 = tmpchi2;
                //StatusCode dump2 = m_dist->distance((const LHCb::Particle *) part,(const LHCb::Vertex *)vd,D_ip,D_chi2);

                if (chi2 < 50) {
                    dummy = 4000;
                    float bdtval = m_Reader->EvaluateMVA("BDT method");
                    for (int i = 0; i < m_nWrite; i++) {
                        if (bdtval > maxbdt[i]) {
                            for (int j = m_nWrite; j > i; j--) {
                                SelParts[j] = SelParts[j - 1];
                                maxbdt[j] = maxbdt[j - 1];
                            }
                            SelParts[i] = part;
                            maxbdt[i] = bdtval;
                            break;
                        }
                    }
                }
            }
        } // end particles loop
    } //end particle types loop
    std::stringstream sout;
    for (int i = 0; i < m_nWrite; i++) {
        std::string name = "";
        if (i > 0) {
            sout << i + 1;
            name = sout.str();
        }
        writeParticle(SelParts[i], maxbdt[i], name, tuple, prefix, mother);
        sout.str("");
    }

    tuple->column(prefix + "_ISOLATION_DstWindowBDT", DstBDT1);
    tuple->column(prefix + "_ISOLATION_DstWindowDELTAM", DstM1);
    tuple->column(prefix + "_ISOLATION_DstWindowBDT2", DstBDT2);
    tuple->column(prefix + "_ISOLATION_DstWindowDELTAM2", DstM2);
    tuple->column(prefix + "_ISOLATION_DstWindowInWindow", nWindow);

    return StatusCode(test);
}

//=========================================================================
//
//=========================================================================
const Vertex* TupleToolApplyIsolationVetoDst::originVertex(const Particle* top, const Particle* P) const
{
    if (top == P || P->isBasicParticle())
        return 0;

    const SmartRefVector<LHCb::Particle>& dau = top->daughters();
    if (dau.empty()) {
        // if (msgLevel(MSG::DEBUG)) debug() << " Particle has no daughters! "  << endmsg;
        return 0;
    }

    SmartRefVector<LHCb::Particle>::const_iterator it;
    for (it = dau.begin(); dau.end() != it; ++it) {
        if (P == *it) { // I found the daughter
            return top->endVertex();
        }
    }

    // vertex not yet found, get deeper in the decay:
    for (it = dau.begin(); dau.end() != it; ++it) {
        if (P != *it && !(*it)->isBasicParticle()) {
            const Vertex* vv = originVertex(*it, P);
            if (vv)
                return vv;
        }
    }
    return 0;
}

//=============================================================================
// Check if the track is already in the decay
//=============================================================================
bool TupleToolApplyIsolationVetoDst::isTrackInDecay(const LHCb::Track* track, std::vector<const LHCb::Track*> daughters)
{
    bool isInDecay = false;
    //loop over daughters
    for (std::vector<const LHCb::Track*>::iterator it = daughters.begin(); it != daughters.end(); ++it) {
        const LHCb::Track* itrack = (*it);
        if (itrack) {
            if (itrack == track) {
                if (msgLevel(MSG::DEBUG))
                    debug() << "Track is in decay, skipping it" << endmsg;
                isInDecay = true;
            }
        }
    } //end daughter loop

    return isInDecay;
}

//=============================================================================
// MINIPCHI2 for a track
//=============================================================================
double TupleToolApplyIsolationVetoDst::getminipchi(const LHCb::Particle* track)
{

    double minchi2 = -1;
    const RecVertex::Range PV = m_dva->primaryVertices();
    if (!PV.empty()) {
        for (RecVertex::Range::const_iterator pv = PV.begin(); pv != PV.end(); ++pv) {
            double ip, chi2;
            StatusCode test2 = m_dist->distance((const LHCb::Particle*)track, *pv, ip, chi2);
            if ((chi2 < minchi2) || (minchi2 < 0.)) {
                LHCb::RecVertex newPV(**pv);
                StatusCode scfit = m_pvReFitter->remove(track, &newPV);
                LHCb::RecVertex* newPVPtr = (LHCb::RecVertex*)&newPV;
                test2 = m_dist->distance((LHCb::Particle*)track, (LHCb::VertexBase*)newPVPtr, ip, chi2);
                minchi2 = chi2;
            }
        }
    }

    return minchi2;
}

double TupleToolApplyIsolationVetoDst::getfdchi2(const LHCb::Track* track, LHCb::Vertex Vtx)
{

    double minchi2 = -1;
    double fdchi2 = -1;
    double fd;
    const RecVertex::Range PV = m_dva->primaryVertices();
    if (!PV.empty()) {
        for (RecVertex::Range::const_iterator pv = PV.begin(); pv != PV.end(); ++pv) {
            double ip, chi2;
            StatusCode test2 = m_dist->distance((const LHCb::Track*)track, *pv, ip, chi2);
            if ((chi2 < minchi2) || (minchi2 < 0.)) {
                minchi2 = chi2;
                StatusCode test2 = m_dist->distance(*pv, &Vtx, fd, fdchi2);
            }
        }
    }

    return fdchi2;
}

//=============================================================================
// Opening angle for a track and particle
//=============================================================================
double TupleToolApplyIsolationVetoDst::getopening(const LHCb::Track* track, const LHCb::Particle* P)
{
    Gaudi::XYZVector A = P->momentum().Vect();
    Gaudi::XYZVector B = track->momentum();
    double cosopening = A.Dot(B) / std::sqrt(A.Mag2() * B.Mag2());
    return cosopening;
}

//=============================================================================
// Store information in tuple for selected particles
//=============================================================================
void TupleToolApplyIsolationVetoDst::writeParticle(const LHCb::Particle* P, double bdt, std::string name, Tuples::Tuple& tuple, std::string prefix, const LHCb::Particle* Mother)
{
    if (bdt > -1) {
        float pe = P->momentum().E();
        float px = P->momentum().Px();
        float py = P->momentum().Py();
        float pz = P->momentum().Pz();
        float pidk = P->proto()->info(LHCb::ProtoParticle::CombDLLk, -1000);
        float pidp = P->proto()->info(LHCb::ProtoParticle::CombDLLp, -1000);
        float nnp = P->proto()->info(LHCb::ProtoParticle::ProbNNp, -1);
        float nnk = P->proto()->info(LHCb::ProtoParticle::ProbNNk, -1);
        float nnpi = P->proto()->info(LHCb::ProtoParticle::ProbNNpi, -1);
        float nng = P->proto()->info(LHCb::ProtoParticle::ProbNNghost, -1);
        float charge = 0;
        if (P->proto()->track()->type() > 2) {
            charge = P->proto()->track()->charge();
        }
        float type = P->proto()->track()->type();
        const MuonPID* muonPID = P->proto()->muonPID();
        float ismuon = muonPID ? muonPID->IsMuon() : false;
        tuple->column(prefix + "_ISOLATION_BDT" + name + m_outputSuffix, bdt);
        tuple->column(prefix + "_ISOLATION_Type" + name + m_outputSuffix, type);
        if (m_verbose) {
            tuple->column(prefix + "_ISOLATION_CHARGE" + name + m_outputSuffix, charge);
            tuple->column(prefix + "_ISOLATION_PE" + name + m_outputSuffix, pe);
            tuple->column(prefix + "_ISOLATION_PX" + name + m_outputSuffix, px);
            tuple->column(prefix + "_ISOLATION_PY" + name + m_outputSuffix, py);
            tuple->column(prefix + "_ISOLATION_PZ" + name + m_outputSuffix, pz);
            tuple->column(prefix + "_ISOLATION_PIDK" + name + m_outputSuffix, pidk);
            tuple->column(prefix + "_ISOLATION_PIDp" + name + m_outputSuffix, pidp);
            tuple->column(prefix + "_ISOLATION_NNk" + name + m_outputSuffix, nnk);
            tuple->column(prefix + "_ISOLATION_NNpi" + name + m_outputSuffix, nnpi);
            tuple->column(prefix + "_ISOLATION_NNp" + name + m_outputSuffix, nnp);
            tuple->column(prefix + "_ISOLATION_IsMuon" + name + m_outputSuffix, ismuon);
            tuple->column(prefix + "_ISOLATION_NNghost" + name + m_outputSuffix, nng);
        }
        if (m_trueID) {
            const LHCb::MCParticle* mcp(NULL);
            const LHCb::MCParticle* mother_mcp(NULL);
            bool foundMCP = false;
            bool foundMCP_mother = false;
            if (P) {
                //assignedPid = P->particleID().pid();
                if (msgLevel(MSG::VERBOSE))
                    verbose() << "Getting related MCP to " << P << endmsg;
                for (std::vector<IParticle2MCAssociator*>::const_iterator iMCAss = m_p2mcAssocs.begin();
                     iMCAss != m_p2mcAssocs.end(); ++iMCAss) {
                    if (!foundMCP)
                        mcp = (*iMCAss)->relatedMCP(P);
                    if (mcp) {
                        foundMCP = true;
                    }
                    if (!foundMCP_mother)
                        mother_mcp = (*iMCAss)->relatedMCP(Mother);
                    if (mother_mcp) {
                        foundMCP_mother = true;
                    }
                }
                if (msgLevel(MSG::VERBOSE))
                    verbose() << "Got mcp " << mcp << endmsg;
            }
            float pid = -1;
            if (foundMCP)
                pid = mcp->particleID().abspid();
            tuple->column(prefix + "_ISOLATION_TRUEPID" + name + m_outputSuffix, pid);
            float motherID = -1;
            if (foundMCP)
                if (mcp->mother())
                    motherID = mcp->mother()->particleID().abspid();
            tuple->column(prefix + "_ISOLATION_MOTHERID" + name + m_outputSuffix, motherID);
            float sameMother = 0;
            if (foundMCP) {
                if (mother_mcp) {
                    const LHCb::MCParticle* head = mcp;
                    for (;;) {
                        head = head->mother();
                        if (!head) {
                            break;
                        }
                        int motherid = head->particleID().abspid();
                        if ((head == mother_mcp) || (head->key() == mother_mcp->key()))
                            sameMother = 1;
                    }
                }
            }

            tuple->column(prefix + "_ISOLATION_SAMEMOTHER" + name + m_outputSuffix, sameMother);
        }

    } else {
        tuple->column(prefix + "_ISOLATION_BDT" + name + m_outputSuffix, bdt);
        tuple->column(prefix + "_ISOLATION_Type" + name + m_outputSuffix, (float)0.);
        if (m_verbose) {
            tuple->column(prefix + "_ISOLATION_CHARGE" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_PE" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_PX" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_PY" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_PZ" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_PIDK" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_PIDp" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_NNk" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_NNpi" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_NNp" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_IsMuon" + name + m_outputSuffix, (float)0.);
            tuple->column(prefix + "_ISOLATION_NNghost" + name + m_outputSuffix, (float)0.);
        }
        if (m_trueID) {
            tuple->column(prefix + "_ISOLATION_TRUEPID" + name + m_outputSuffix, (float)-1.);
            tuple->column(prefix + "_ISOLATION_MOTHERID" + name + m_outputSuffix, (float)-1.);
            tuple->column(prefix + "_ISOLATION_SAMEMOTHER" + name + m_outputSuffix, (float)-1);
        }
    }
}
