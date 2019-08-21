// $Id: TupleToolTauMuDiscrVars.cpp,v 1.17 2010-05-12 20:01:40 jpalac Exp $
// Include files

// local
#include "TupleToolTauMuDiscrVars.h"

#include "Kernel/IDVAlgorithm.h"
#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>
#include "Kernel/IPVReFitter.h"

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/Particle.h"

#include "LoKi/Math.h"

using namespace LHCb;

//-----------------------------------------------------------------------------
// Implementation file for class : GeometryTupleTool
//
// 2007-11-07 : Jeremie Borel
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_TOOL_FACTORY( TupleToolTauMuDiscrVars )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
  TupleToolTauMuDiscrVars::TupleToolTauMuDiscrVars( const std::string& type,
                                        const std::string& name,
                                        const IInterface* parent )
    :
    TupleToolBase ( type, name , parent )
    , m_dist(0)
    , m_dva(0)
    , m_pvReFitterName( "LoKi::PVReFitter:PUBLIC" )
{
  declareInterface<IParticleTupleTool>(this);
  declareProperty("RefitPVs",m_refitPVs=false,
                  "Refit PVs when doing next best PV checks");
  declareProperty("PVReFitter", m_pvReFitterName,
                  "PV refitter algorithm name (':PUBLIC' at end of algo name makes sure a public instance is used)" );
  declareProperty("FillMultiPV",m_fillMultiPV=false,
                  "Fill Multi PV arrays");

  //declareProperty("FillMother",m_fillMother=true,
  //                "Turn false if the mother is expected to be NULL, will not fill mother PV info");
  // replaced by Verbose

}

//=============================================================================

StatusCode TupleToolTauMuDiscrVars::initialize()
{
  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc(), this ) ;
  if (!m_dva) return Error("Couldn't get parent DVAlgorithm");

  m_dist = m_dva->distanceCalculator();
  if ( !m_dist )
  {
    return Error("Unable to retrieve the IDistanceCalculator tool");
  }

  m_pvReFitter = tool<IPVReFitter>( m_pvReFitterName, this );
  if ( !m_pvReFitter )
  {
    return Error( "Unable to retrieve IPVReFitter instance" );
  }

  return sc;
}

//=============================================================================
StatusCode TupleToolTauMuDiscrVars::fill( const Particle* mother
                                    , const Particle* P
                                    , const std::string& head
                                    , Tuples::Tuple& tuple )
{

  bool test = true;

  if( mother && P ){

    // Primary and secondary vertices.
    const VertexBase* B_PV = m_dva->bestVertex(P);
    const VertexBase* B_SV = P->endVertex();

    // Bd flight direction.
    const LoKi::ThreeVector B_direction = B_SV->position() - B_PV->position();
    const double norm = ::sqrt(B_direction.Mag2());
    const double B_cosangle_X = B_direction.X()/norm;
    const double B_cosangle_Y = B_direction.Y()/norm;
    const double B_cosangle_Z = B_direction.Z()/norm;
    const double B_angle_Z = ::acos(B_cosangle_Z);

    // Bd estimated momentum.
    const double B_M = 5279.61;
    const double observed_B_M = P->measuredMass();
    const double observed_B_PZ = P->momentum().Z();
    const double B_P = B_M/observed_B_M*observed_B_PZ*::sqrt(1.+::tan(B_angle_Z)*::tan(B_angle_Z));
    const LoKi::LorentzVector B_PXYZE (B_P*B_cosangle_X,B_P*B_cosangle_Y,B_P*B_cosangle_Z,::sqrt(B_P*B_P+B_M*B_M));

    // Bd daughters momenta.
    const LoKi::LorentzVector D_PXYZE = P->daughters()[0]->momentum();
    const LoKi::LorentzVector mu_PXYZE = P->daughters()[1]->momentum();

    // Bd boost.
    const double B_beta = B_P/::sqrt(B_P*B_P+B_M*B_M);
    const double bx = -B_beta*B_cosangle_X;
    const double by = -B_beta*B_cosangle_Y;
    const double bz = -B_beta*B_cosangle_Z;
    const double b2 = B_beta*B_beta;
    const double gamma = 1./::sqrt(1.-b2);
    const double bp = bx*mu_PXYZE.X()+by*mu_PXYZE.Y()+bz*mu_PXYZE.Z();
    const double gamma2 = b2>0?(gamma-1.)/b2:0.;

    // Muon momentum in the Bd rest frame.
    const double mu_PXYZE_Brestframe_X = mu_PXYZE.X() + gamma2*bp*bx + gamma*bx*mu_PXYZE.E();
    const double mu_PXYZE_Brestframe_Y = mu_PXYZE.Y() + gamma2*bp*by + gamma*by*mu_PXYZE.E();
    const double mu_PXYZE_Brestframe_Z = mu_PXYZE.Z() + gamma2*bp*bz + gamma*bz*mu_PXYZE.E();
    const double mu_PXYZE_Brestframe_E = gamma*(mu_PXYZE.E() + bp);
    const LoKi::LorentzVector mu_PXYZE_Brestframe (mu_PXYZE_Brestframe_X,mu_PXYZE_Brestframe_Y,mu_PXYZE_Brestframe_Z,mu_PXYZE_Brestframe_E);

    // Discriminant variables.
    const double El = mu_PXYZE_Brestframe.E();
    const LoKi::LorentzVector Pcomb1 = B_PXYZE-D_PXYZE-mu_PXYZE;
    const double Mmiss2 = Pcomb1.E()*Pcomb1.E()-Pcomb1.X()*Pcomb1.X()-Pcomb1.Y()*Pcomb1.Y()-Pcomb1.Z()*Pcomb1.Z();
    const LoKi::LorentzVector Pcomb2 = B_PXYZE-D_PXYZE;
    const double q2 = Pcomb2.E()*Pcomb2.E()-Pcomb2.X()*Pcomb2.X()-Pcomb2.Y()*Pcomb2.Y()-Pcomb2.Z()*Pcomb2.Z();

    // Filling of the branches.
    test &= tuple->column( head+"_FlightDir_Zangle" , B_angle_Z );
    test &= tuple->column( head+"_Estimated_P" , B_P );
    test &= tuple->column( head+"_Estimated_PX" , B_PXYZE.X() );
    test &= tuple->column( head+"_Estimated_PY" , B_PXYZE.Y() );
    test &= tuple->column( head+"_Estimated_PZ" , B_PXYZE.Z() );
    test &= tuple->column( head+"_Boost_Beta" , B_beta );
    test &= tuple->column( "FitVar_El" , El );
    test &= tuple->column( "FitVar_Mmiss2" , Mmiss2 );
    test &= tuple->column( "FitVar_q2" , q2 );

  }

  else 
  {
    return StatusCode::FAILURE;
  }

  return StatusCode(test);

}

//=========================================================================
//  Fill Everything for this vertex for related PV
//=========================================================================
StatusCode TupleToolTauMuDiscrVars::fillVertexFull(const LHCb::VertexBase* vtx,
                                             const LHCb::Particle* P,
                                             const std::string& prefix,
                                             const std::string& vtx_name,
                                             Tuples::Tuple& tuple) const
{
  if ( !vtx ) ++counter("Can't retrieve the " +vtx_name+ " vertex for " + prefix );
  StatusCode sc = fillVertex(vtx,prefix+vtx_name,tuple);
  if ( sc.isFailure() )
  {
    return Warning("Could not fill Endvertex "+prefix,sc,1);
  }
  sc = fillBPV(vtx,P,prefix,tuple,vtx_name);
  if ( sc.isFailure() )
  {
    return Warning("Could not fillBPV "+prefix,sc,1);
  }
  if( !P->isBasicParticle() ) sc = fillFlight(vtx,P,prefix,tuple,vtx_name);
  if ( sc.isFailure() ) Warning("Error in fillVertexFull "+prefix,sc,1).ignore();
  return sc;
}

//=========================================================================
//  Fill PV for related PV
//=========================================================================
StatusCode TupleToolTauMuDiscrVars::fillBPV( const VertexBase* primVtx
                                       , const Particle* P
                                       , const std::string& prefix
                                       , Tuples::Tuple& tuple
                                       , const std::string& trail) const {
  bool test = true ;

  double ip=0, chi2=0;
  if ( !primVtx )
  {
    ++counter("No BPV for "+prefix);
    test &= tuple->column( prefix + "_IP"+trail, -999. );
    test &= tuple->column( prefix + "_IPCHI2"+trail, -999. );
  }
  else
  {
    test &= m_dist->distance ( P, primVtx, ip, chi2 );
    if ( !test )
    {
      ip   = -1;
      chi2 = -1;
    }
    test &= tuple->column( prefix + "_IP"+trail, ip );
    test &= tuple->column( prefix + "_IPCHI2"+trail, chi2 );
  }

  if (!test) Warning("Error in fillBPV "+prefix,StatusCode(test),1).ignore();
  return StatusCode(test) ;
}
//=========================================================================
//  Fill PV for all PV
//=========================================================================
StatusCode TupleToolTauMuDiscrVars::fillMinIP( const Particle* P,
                                         const std::string& prefix,
                                         Tuples::Tuple& tuple ) const
{
  bool test = true ;
  // minimum IP
  double ipmin = -1;
  double minchi2 = -1 ;

  double ipminnextbest = -1;
  double minchi2nextbest = -1;
  if(msgLevel(MSG::VERBOSE)) verbose() << "Looking for Min IP"  << endmsg  ;
  const RecVertex::Range PV = m_dva->primaryVertices();
  if(msgLevel(MSG::VERBOSE)) verbose() << "PV size: "  << PV.size() << endmsg  ;
  std::vector<double> ips, ipchi2s, diras;
  if ( !PV.empty() )
  {
    if(msgLevel(MSG::VERBOSE)) verbose() << "Filling IP " << prefix + "_MINIP : "
                                         << P << " PVs:" << PV.size() << endmsg ;

    for ( RecVertex::Range::const_iterator pv = PV.begin() ; pv!=PV.end() ; ++pv)
    {
      RecVertex newPV(**pv);
      if (m_refitPVs)
      {

        StatusCode scfit = m_pvReFitter->remove(P, &newPV);
        if(!scfit) { Warning("ReFitter fails!",StatusCode::SUCCESS,10).ignore(); continue; }
      }

      double ip, chi2;
      //StatusCode test2 = m_dist->distance ( P, *pv, ip, chi2 );

      LHCb::VertexBase* newPVPtr = (LHCb::VertexBase*)&newPV;
      StatusCode test2 = m_dist->distance ( P, newPVPtr, ip, chi2 );
      ips.push_back(ip);
      ipchi2s.push_back(chi2);
      if (P->endVertex()) diras.push_back(dira(newPVPtr,P));
      if( test2 && isVerbose() )
      {
        if ((ip<ipmin) || (ipmin<0.))
        {
          ipminnextbest = ipmin;
          ipmin = ip ;
        }
        else
        {
          if((ip < ipminnextbest) || (ipminnextbest < 0))
          {
            ipminnextbest = ip;
          }
        }

        if ((chi2<minchi2) || (minchi2<0.))
        {
          minchi2nextbest = minchi2;
          minchi2 = chi2 ;
        }
        else
        {
          if((chi2 < minchi2nextbest) || (minchi2nextbest < 0))
          {
            minchi2nextbest = chi2;
          }
        }
      }
    }
  }
  if (isVerbose()){
    if ( msgLevel(MSG::VERBOSE) )
    {
      verbose() << "Filling IP " << prefix + "_MINIP " << ipmin << " at " << minchi2 << endmsg  ;
      verbose() << "Filling IP next best " << prefix + "_MINIP " << ipminnextbest << " at "
                << minchi2nextbest << endmsg  ;
    }
    test &= tuple->column( prefix + "_MINIP", ipmin );
    test &= tuple->column( prefix + "_MINIPCHI2", minchi2 );
    
    test &= tuple->column( prefix + "_MINIPNEXTBEST", ipminnextbest );
    test &= tuple->column( prefix + "_MINIPCHI2NEXTBEST", minchi2nextbest );
  }
  if (m_fillMultiPV){
    test &= tuple->farray( prefix + "_AllIP", ips, "nPV", m_maxPV );
    test &= tuple->farray( prefix + "_AllIPchi2", ipchi2s, "nPV", m_maxPV );
    if (!diras.empty()) test &= tuple->farray( prefix + "_AllDIRA", diras, "nPV", m_maxPV );
    // --------------------------------------------------
  }
  
  if(msgLevel(MSG::VERBOSE))
    verbose() << "Return from fillMinIP: " << prefix  << " " << test << endmsg;
  if (!test)  Warning("Error in fillMinIP",StatusCode::FAILURE,1);
  return StatusCode(test) ;
}
//=========================================================================
// fill vertex stuff
//=========================================================================
StatusCode TupleToolTauMuDiscrVars::fillVertex( const LHCb::VertexBase* vtx,
                                          const std::string& vtx_name,
                                          Tuples::Tuple& tuple ) const
{
  bool test = true ;

  // decay vertex information:
  if ( !vtx )
  {
    Gaudi::XYZPoint pt(-999.,-999.,-999.) ; // arbitrary point
    test &= tuple->column(  vtx_name+"_", pt );
    test &= tuple->column(  vtx_name + "_XERR", -999. );
    test &= tuple->column(  vtx_name + "_YERR", -999. );
    test &= tuple->column(  vtx_name + "_ZERR", -999. );
    test &= tuple->column(  vtx_name + "_CHI2", -999. );
    test &= tuple->column(  vtx_name + "_NDOF", -1 );
    test &= tuple->matrix(  vtx_name + "_COV_", Gaudi::SymMatrix3x3()  );
  }
  else
  {
    test &= tuple->column( vtx_name+"_", vtx->position() );
    const Gaudi::SymMatrix3x3 & m = vtx->covMatrix ();
    test &= tuple->column(  vtx_name + "_XERR", std::sqrt( m(0,0) ) );
    test &= tuple->column(  vtx_name + "_YERR", std::sqrt( m(1,1) ) );
    test &= tuple->column(  vtx_name + "_ZERR", std::sqrt( m(2,2) ) );
    test &= tuple->column(  vtx_name + "_CHI2", vtx->chi2() );
    test &= tuple->column(  vtx_name + "_NDOF", vtx->nDoF() );
    test &= tuple->matrix(  vtx_name + "_COV_", m );
  }

  // --------------------------------------------------
  if (!test) Warning("Error in fillVertex "+vtx_name,StatusCode(test),1).ignore();
  return StatusCode(test) ;

}
//=========================================================================
// fill flight distance, angle...
//=========================================================================
StatusCode TupleToolTauMuDiscrVars::fillFlight( const VertexBase* oriVtx,
                                          const Particle* P,
                                          const std::string& prefix,
                                          Tuples::Tuple& tuple,
                                          const std::string& trail ) const
{
  bool test = true ;
  // --------------------------------------------------
  if ( !oriVtx )
  {
    test &= tuple->column( prefix + "_FD"+trail, -999. );
    test &= tuple->column( prefix + "_FDCHI2"+trail, -999. );
    test &= tuple->column( prefix + "_DIRA"+trail, -999.);
  }
  else
  {

    // flight distance
    double dist = 0;
    double chi2 = 0 ;
    StatusCode sc = m_dist->distance( oriVtx, P->endVertex(), dist, chi2 );
    if ( sc.isFailure() ) return sc ;

    test &= tuple->column( prefix + "_FD"+trail, dist );
    test &= tuple->column( prefix + "_FDCHI2"+trail, chi2 );
    // --------------------------------------------------
    // cosine of (flight distance) dot (momentum):
    // find the origin vertex. Either the primary or the origin in the
    // decay
    test &= tuple->column( prefix + "_DIRA"+trail, dira(oriVtx,P) );
  }

  if (!test) Warning("Error in fillFlight "+prefix,StatusCode(test),1).ignore();
  return StatusCode(test);
}
// =====================================================
// find origin vertex in the decay chain
// =====================================================
const VertexBase* TupleToolTauMuDiscrVars::originVertex( const Particle* top,
                                                   const Particle* P ) const
{
  //this used to pass back zero if P was a basic particle.
  //I don't think that's necessary. R Lambert 2009-08-14
  if( top == P || top->isBasicParticle() ) return 0;

  const SmartRefVector<LHCb::Particle>& dau = top->daughters ();
  if ( dau.empty() ) return 0;

  for ( SmartRefVector<LHCb::Particle>::const_iterator it = dau.begin();
        dau.end() != it; ++it )
  {
    if( P == *it )
    { // I found the daughter
      if(msgLevel(MSG::VERBOSE))
        verbose() << "It's a daughter, retrning mother's endvertex : "
                  << endmsg;
      return top->endVertex();
    }
  }

  // vertex not yet found, get deeper in the decay:
  for( SmartRefVector<LHCb::Particle>::const_iterator it = dau.begin();
       dau.end() != it; ++it )
  {
    if ( P != *it && !(*it)->isBasicParticle() )
    {
      const VertexBase* vv = originVertex( *it, P );
      if ( msgLevel(MSG::VERBOSE) ) verbose() << "Went up : " << vv  << endmsg  ;
      if( vv ) return vv;
    }
  }
  return 0;
}
