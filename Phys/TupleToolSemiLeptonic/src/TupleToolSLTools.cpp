// local
#include "TupleToolSLTools.h"

// from Phys
#include "Kernel/GetIDVAlgorithm.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IDistanceCalculator.h"
#include "Kernel/ILifetimeFitter.h"
#include "Kernel/IPVReFitter.h"
#include "Kernel/IVertexFit.h"

// from Gaudi
#include "GaudiAlg/Tuple.h"

// from LHCb
#include "Event/Particle.h"

// from ROOT
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TVector3.h"

using namespace LHCb;

DECLARE_COMPONENT( TupleToolSLTools )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolSLTools::TupleToolSLTools( const std::string& type,
                                    const std::string& name,
                                    const IInterface*  parent )
    : TupleToolBase( type, name, parent ),
      m_dva( 0 ),
      m_dist( 0 ),
      m_pVertexFit( 0 ) {
  declareInterface<IParticleTupleTool>( this );
  declareProperty( "VertexFitter", m_typeVertexFit = "LoKi::VertexFitter" );
  declareProperty( "Bmass", m_Bmass = 5619.5 );
  declareProperty( "VertexCov", m_vcov = false );
  declareProperty( "MomCov", m_momcov = false );
}

//=============================================================================
StatusCode TupleToolSLTools::initialize() {
  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc(), this );
  if ( !m_dva )
    return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE );

  m_dist = m_dva->distanceCalculator();
  if ( !m_dist ) {
    Error( "Unable to retrieve the IDistanceCalculator tool" );
    return StatusCode::FAILURE;
  }

  // m_pVertexFit= m_dva->vertexFitter();
  m_pVertexFit = tool<IVertexFit>( m_typeVertexFit, this );
  // m_pVertexFit  = tool<IVertexFit>("OfflineVertexFitter",this);
  if ( !m_pVertexFit ) {
    Error( "Unable to retrieve the IVertexFit tool" );
    return StatusCode::FAILURE;
  }

  return sc;
}

//=============================================================================
StatusCode TupleToolSLTools::fill( const Particle* mother, const Particle* P,
                                   const std::string& head,
                                   Tuples::Tuple&     tuple ) {
  const std::string prefix = fullName( head );
  Assert( P && mother && m_dist,
          "This should not happen, you are inside TupleToolSLTools.cpp :(" );

  bool test = true;

  // find the origin vertex. Either the primary or the origin in the
  // decay
  const LHCb::Particle*       mu_part( nullptr );
  LHCb::Particle::ConstVector source;
  LHCb::Particle::ConstVector target;

  double                      mcorr, q2_one, q2_two, mcorrerr, mcorrfullerr;
  std::vector<double>         mcorr_errors;
  std::vector<TLorentzVector> nu_slns;
  const LHCb::Vertex*         pmu_vert;
  const LHCb::VertexBase*     PV;

  if ( P->isBasicParticle() ) {
    source.push_back( mother );
  } else {
    source.push_back( P );
  }
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
            if ( ( *itmp )->particleID().abspid() == 13 ) {
              mu_part = ( *itmp );
            }
          }
        }
      }  // if endVertex
    }    // isource
    source = target;
  } while ( !target.empty() );

  // vetex selected proton and muon
  pmu_vert = P->endVertex();
  PV       = m_dva->bestVertex( P );

  Gaudi::LorentzVector LV_P, LV_mu;
  LV_P = P->momentum();
  TLorentzVector TLV_P( LV_P.px(), LV_P.py(), LV_P.pz(), LV_P.E() );
  // muon LV
  LV_mu = mu_part->momentum();
  TLorentzVector TLV_mu( LV_mu.px(), LV_mu.py(), LV_mu.pz(), LV_mu.E() );

  TVector3 TV3_PV( PV->position().X(), PV->position().Y(),
                   PV->position().Z() );
  TVector3 TV3_SV( pmu_vert->position().X(), pmu_vert->position().Y(),
                   pmu_vert->position().Z() );
  // double dist_test = (TV3_SV-TV3_PV).Mag();
  TVector3 Mdirn = ( TV3_SV - TV3_PV ).Unit();
  mcorr          = Mcorr( TLV_P, Mdirn );
  mcorr_errors =
      Mcorr_errors( TV3_SV, TV3_PV, TLV_P, P->covMatrix(), PV->covMatrix() );
  std::vector<double> q2_err = q2_errors( TV3_SV, TV3_PV, TLV_P, TLV_mu,
                                          P->covMatrix(), PV->covMatrix() );
  nu_slns                    = recoNu( TLV_P, Mdirn, m_Bmass );
  fillVertex( pmu_vert, prefix + "_SV", tuple );
  fillVertex( PV, prefix + "_PV", tuple );
  double         dist    = 0;
  double         chi2    = 0;
  StatusCode     sc      = m_dist->distance( PV, P->endVertex(), dist, chi2 );
  double         theta   = Mdirn.Theta();
  double         InvSinF = 1.0 / sin( theta );
  double         Preg_B  = 1000 * ( InvSinF * 1.8004 + dist * 4.95956 );
  TLorentzVector regBfourvect(
      Preg_B * Mdirn.X(), Preg_B * Mdirn.Y(), Preg_B * Mdirn.Z(),
      sqrt( Preg_B * Mdirn.Z() * Preg_B * Mdirn.Z() + m_Bmass * m_Bmass ) );
  double Preg = ( regBfourvect - TLV_P ).P();

  if ( m_momcov ) {
    const Gaudi::SymMatrix4x4& mom_cov    = P->momCovMatrix();
    const Gaudi::Matrix4x3&    posmom_cov = P->posMomCovMatrix();
    test &= tuple->matrix( prefix + "_MOM_COV_", mom_cov );
    test &= tuple->matrix( prefix + "_POSMOM_COV_", posmom_cov );
  }

  test &= tuple->column( prefix + "_MCORR", mcorr );
  if ( q2_err.size() == 2 ) {
    test &= tuple->column( prefix + "_Q2ERR_LOW", q2_err[0] );
    test &= tuple->column( prefix + "_Q2ERR_HIGH", q2_err[1] );
  } else {
    test &= tuple->column( prefix + "_Q2ERR_LOW", -1000. );
    test &= tuple->column( prefix + "_Q2ERR_HIGH", -1000. );
  }
  if ( mcorr_errors.size() == 2 ) {
    mcorrerr     = mcorr_errors[0];
    mcorrfullerr = mcorr_errors[1];
    test &= tuple->column( prefix + "_MCORRERR", mcorrerr );
    test &= tuple->column( prefix + "_MCORRFULLERR", mcorrfullerr );
  } else {
    test &= tuple->column( prefix + "_MCORRERR", -1000. );
    test &= tuple->column( prefix + "_MCORRFULLERR", -1000. );
  }

  if ( nu_slns.size() == 2 ) {
    TLorentzVector munu_one = nu_slns[0] + TLV_mu;
    TLorentzVector munu_two = nu_slns[1] + TLV_mu;
    TLorentzVector Lb_one   = nu_slns[0] + TLV_P;
    TLorentzVector Lb_two   = nu_slns[1] + TLV_P;

    TVector3 boostvec_one = -munu_one.BoostVector();
    TVector3 boostvec_two = -munu_two.BoostVector();

    TLorentzVector TLV_mu_one( TLV_mu );
    TLorentzVector TLV_mu_two( TLV_mu );

    Lb_one.Boost( boostvec_one );
    TLV_mu_one.Boost( boostvec_one );
    Lb_two.Boost( boostvec_two );
    TLV_mu_two.Boost( boostvec_two );

    double costhetal_one = -cos( Lb_one.Angle( TLV_mu_one.Vect() ) );
    double costhetal_two = -cos( Lb_two.Angle( TLV_mu_two.Vect() ) );

    q2_one = ( nu_slns[0] + TLV_mu ).M2();
    q2_two = ( nu_slns[1] + TLV_mu ).M2();

    if ( fabs( nu_slns[0].P() - Preg ) < fabs( nu_slns[1].P() - Preg ) ) {
      test &= tuple->column( prefix + "_Q2BEST", q2_one );
      test &= tuple->column( prefix + "_COSTHETALBEST", costhetal_one );
      test &= tuple->column( prefix + "_Q2BESTERR", q2_err[0] );
    } else {
      test &= tuple->column( prefix + "_Q2BEST", q2_two );
      test &= tuple->column( prefix + "_COSTHETALBEST", costhetal_two );
      test &= tuple->column( prefix + "_Q2BESTERR", q2_err[1] );
    }
    test &= tuple->column( prefix + "_Q2SOL1", q2_one );
    test &= tuple->column( prefix + "_Q2SOL2", q2_two );
    test &= tuple->column( prefix + "_COSTHETALSOL1", costhetal_one );
    test &= tuple->column( prefix + "_COSTHETALSOL2", costhetal_two );
    test &= tuple->column( prefix + "_PNUSOL1", nu_slns[0].P() );
    test &= tuple->column( prefix + "_PNUSOL2", nu_slns[1].P() );
    test &= tuple->column( prefix + "_PREG", Preg );
  } else {
    test &= tuple->column( prefix + "_Q2SOL1", -1000. );
    test &= tuple->column( prefix + "_Q2SOL2", -1000. );
    test &= tuple->column( prefix + "_COSTHETALSOL1", -1000. );
    test &= tuple->column( prefix + "_COSTHETALSOL2", -1000. );
    test &= tuple->column( prefix + "_PNUSOL1", -1000. );
    test &= tuple->column( prefix + "_PNUSOL2", -1000. );
    test &= tuple->column( prefix + "_PREG", -1000. );
    test &= tuple->column( prefix + "_Q2BEST", -1000. );
    test &= tuple->column( prefix + "_COSTHETALBEST", -1000. );
    test &= tuple->column( prefix + "_Q2BESTERR", -1000. );
  }
  debug() << "PAST FILL " << endmsg;

  return StatusCode( test );
}

//=========================================================================
StatusCode TupleToolSLTools::fillVertex( const LHCb::VertexBase* vtx,
                                         const std::string&      vtx_name,
                                         Tuples::Tuple& tuple ) const {
  bool test = true;

  // decay vertex information:
  if ( !vtx ) {
    Gaudi::XYZPoint pt( -999., -999., -999. );  // arbitrary point
    test &= tuple->column( vtx_name + "_", pt );
    test &= tuple->column( vtx_name + "_XERR", -999. );
    test &= tuple->column( vtx_name + "_YERR", -999. );
    test &= tuple->column( vtx_name + "_ZERR", -999. );
    test &= tuple->column( vtx_name + "_CHI2", -999. );
    test &= tuple->column( vtx_name + "_NDOF", -1 );
    test &= tuple->matrix( vtx_name + "_COV_", Gaudi::SymMatrix3x3() );
  } else {
    test &= tuple->column( vtx_name + "_", vtx->position() );
    test &= tuple->column( vtx_name + "_CHI2", vtx->chi2() );
    test &= tuple->column( vtx_name + "_NDOF", vtx->nDoF() );
    if ( m_vcov ) {
      const Gaudi::SymMatrix3x3& m = vtx->covMatrix();
      test &= tuple->matrix( vtx_name + "_COV_", m );
    }
  }
  if ( !test )
    Warning( "Error in fillVertex " + vtx_name, StatusCode( test ), 1 )
        .ignore();
  return StatusCode( test );
}

double TupleToolSLTools::Mcorr( TLorentzVector Y, TVector3 M_dirn ) {
  double P_T = ( Y.Vect() ).Perp( M_dirn );
  return sqrt( Y.M2() + P_T * P_T ) + P_T;
}

std::vector<TLorentzVector> TupleToolSLTools::recoNu( TLorentzVector Y,
                                                      TVector3       M_dirn,
                                                      double         mass ) {
  std::vector<TLorentzVector> p_nu;
  // Get 3 Vector from Y
  TVector3 Y3 = ( Y ).Vect();
  // Construct combined pmu/Kmu 4 vector in new rotated frame
  TLorentzVector Ppmu( Y3.Dot( M_dirn ), ( Y ).Perp( M_dirn ), 0.0,
                       ( Y ).E() );
  // Get component of Y3 perpendicular to *M_dirn
  TVector3 Perp_dirn = ( Y3 - (M_dirn)*Y3.Dot( ( M_dirn ) ) ).Unit();
  // Calculate neutrino energy in mother rest frame
  double E_nurest = ( mass * mass - ( Y ).M2() ) / ( 2 * mass );
  // Calculate pmu " "
  double E_pmurest  = sqrt( E_nurest * E_nurest + ( Y ).M2() );
  double px_rest_sq = ( E_nurest * E_nurest - Ppmu.Py() * Ppmu.Py() );

  // Find magnitude of momentum along mother dirn in rest frame
  // px_rest = sqrt( E_nurest * E_nurest  - Ppmu.Py() * Ppmu.Py() );

  // quadratic coefficients
  double A = ( Y ).E() * ( Y ).E() - Ppmu.Px() * Ppmu.Px();
  double B = -2 * Ppmu.Px() * ( E_nurest * E_pmurest + px_rest_sq );
  double C = -( E_nurest * E_pmurest + px_rest_sq ) *
                 ( E_nurest * E_pmurest + px_rest_sq ) +
             ( Y ).E() * ( Y ).E() * Ppmu.Py() * Ppmu.Py();

  double p1nu;
  double p2nu;
  if ( B * B < 4 * A * C ) {
    p1nu = ( -B ) / ( 2 * A );
    p2nu = p1nu;
  } else {
    // Two neutrino E/p solutions in Lab Frame
    p1nu = ( -B + sqrt( B * B - 4 * A * C ) ) / ( 2 * A );
    p2nu = ( -B - sqrt( B * B - 4 * A * C ) ) / ( 2 * A );
  }

  // reconstruct neutrino 3 vectors and 4 vectors
  TVector3 P_nu_recon_3V1 = (M_dirn)*p1nu + Perp_dirn * -Ppmu.Py();
  TVector3 P_nu_recon_3V2 = (M_dirn)*p2nu + Perp_dirn * -Ppmu.Py();
  p_nu.emplace_back( P_nu_recon_3V1,
                     sqrt( p1nu * p1nu + Ppmu.Py() * Ppmu.Py() ) );
  p_nu.emplace_back( P_nu_recon_3V2,
                     sqrt( p2nu * p2nu + Ppmu.Py() * Ppmu.Py() ) );
  return p_nu;
}

std::vector<double> TupleToolSLTools::Mcorr_errors(
    TVector3 v1, TVector3 v2, TLorentzVector p, Gaudi::SymMatrix7x7 cov1,
    Gaudi::SymMatrix3x3 cov2 ) {
  double x  = v1.Px();
  double y  = v1.Py();
  double z  = v1.Pz();
  double xp = v2.Px();
  double yp = v2.Py();
  double zp = v2.Pz();
  double px = p.Px();
  double py = p.Py();
  double pz = p.Pz();
  double k  = p.E();

  double dMcdpx =
      ( 2 *
            ( 1 - pow( x - xp, 2 ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                       pow( z - zp, 2 ) ) ) *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) -
        ( 2 * ( x - xp ) * ( y - yp ) *
          ( py - ( ( y - yp ) *
                   ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                     ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                       pow( z - zp, 2 ) ) ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
        ( 2 * ( x - xp ) *
          ( pz -
            ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
              ( z - zp ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) /
          ( 2. *
            sqrt(
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) ) +
      ( -2 * px +
        2 * ( 1 - pow( x - xp, 2 ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) -
        ( 2 * ( x - xp ) * ( y - yp ) *
          ( py - ( ( y - yp ) *
                   ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                     ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                       pow( z - zp, 2 ) ) ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
        ( 2 * ( x - xp ) *
          ( pz -
            ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
              ( z - zp ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) /
          ( 2. *
            sqrt(
                pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) - pow( pz, 2 ) +
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) );

  double dMcdpy =
      ( 2 *
            ( 1 - pow( y - yp, 2 ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                       pow( z - zp, 2 ) ) ) *
            ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) -
        ( 2 * ( x - xp ) * ( y - yp ) *
          ( px - ( ( x - xp ) *
                   ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                     ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                       pow( z - zp, 2 ) ) ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
        ( 2 * ( y - yp ) *
          ( pz -
            ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
              ( z - zp ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) /
          ( 2. *
            sqrt(
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) ) +
      ( -2 * py +
        2 * ( 1 - pow( y - yp, 2 ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) -
        ( 2 * ( x - xp ) * ( y - yp ) *
          ( px - ( ( x - xp ) *
                   ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                     ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                       pow( z - zp, 2 ) ) ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
        ( 2 * ( y - yp ) *
          ( pz -
            ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
              ( z - zp ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) /
          ( 2. *
            sqrt(
                pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) - pow( pz, 2 ) +
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) );

  double dMcdpz =
      ( 2 *
            ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                     ( z - zp ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( 1 - pow( z - zp, 2 ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                       pow( z - zp, 2 ) ) ) -
        ( 2 * ( x - xp ) *
          ( px -
            ( ( x - xp ) *
              ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
        ( 2 * ( y - yp ) *
          ( py -
            ( ( y - yp ) *
              ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) /
          ( 2. *
            sqrt(
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) ) +
      ( -2 * pz +
        2 * ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( 1 - pow( z - zp, 2 ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                       pow( z - zp, 2 ) ) ) -
        ( 2 * ( x - xp ) *
          ( px -
            ( ( x - xp ) *
              ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
        ( 2 * ( y - yp ) *
          ( py -
            ( ( y - yp ) *
              ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
          ( z - zp ) ) /
            ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) /
          ( 2. *
            sqrt(
                pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) - pow( pz, 2 ) +
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) );

  double dMcdE =
      k /
      sqrt(
          pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) - pow( pz, 2 ) +
          pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                     pz * ( z - zp ) ) ) /
                        ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                          pow( z - zp, 2 ) ),
               2 ) +
          pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                     pz * ( z - zp ) ) ) /
                        ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                          pow( z - zp, 2 ) ),
               2 ) +
          pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                      ( z - zp ) ) /
                        ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                          pow( z - zp, 2 ) ),
               2 ) );

  double dMcdx = ( 2 *
                       ( ( 2 * pow( x - xp, 2 ) *
                           ( px * ( x - xp ) + py * ( y - yp ) +
                             pz * ( z - zp ) ) ) /
                             pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ),
                                  2 ) -
                         ( px * ( x - xp ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) -
                         ( px * ( x - xp ) + py * ( y - yp ) +
                           pz * ( z - zp ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) *
                       ( px -
                         ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                          pz * ( z - zp ) ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) +
                   2 *
                       ( ( 2 * ( x - xp ) * ( y - yp ) *
                           ( px * ( x - xp ) + py * ( y - yp ) +
                             pz * ( z - zp ) ) ) /
                             pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ),
                                  2 ) -
                         ( px * ( y - yp ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) *
                       ( py -
                         ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                          pz * ( z - zp ) ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) +
                   2 *
                       ( ( 2 * ( x - xp ) *
                           ( px * ( x - xp ) + py * ( y - yp ) +
                             pz * ( z - zp ) ) *
                           ( z - zp ) ) /
                             pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ),
                                  2 ) -
                         ( px * ( z - zp ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) *
                       ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                  pz * ( z - zp ) ) *
                                ( z - zp ) ) /
                                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                    pow( z - zp, 2 ) ) ) ) /
                     ( 2. * sqrt(
                                pow( px - ( ( x - xp ) *
                                            ( px * ( x - xp ) +
                                              py * ( y - yp ) +
                                              pz * ( z - zp ) ) ) /
                                              ( pow( x - xp, 2 ) +
                                                pow( y - yp, 2 ) +
                                                pow( z - zp, 2 ) ),
                                     2 ) +
                                pow( py - ( ( y - yp ) * ( px * ( x - xp ) +
                                                           py * ( y - yp ) + pz * ( z - zp ) ) ) /
                                              ( pow( x - xp, 2 ) +
                                                pow( y - yp, 2 ) +
                                                pow( z - zp, 2 ) ),
                                     2 ) +
                                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                                            ( z - zp ) ) /
                                              ( pow( x - xp, 2 ) +
                                                pow( y - yp, 2 ) +
                                                pow( z - zp, 2 ) ),
                                     2 ) ) ) +
                 ( 2 *
                       ( ( 2 * pow( x - xp, 2 ) *
                           ( px * ( x - xp ) + py * ( y - yp ) +
                             pz * ( z - zp ) ) ) /
                             pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ),
                                  2 ) -
                         ( px * ( x - xp ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) -
                         ( px * ( x - xp ) + py * ( y - yp ) +
                           pz * ( z - zp ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) *
                       ( px -
                         ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                          pz * ( z - zp ) ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) +
                   2 * ( ( 2 * ( x - xp ) * ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) - ( px * ( y - yp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
                       ( py -
                         ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                          pz * ( z - zp ) ) ) /
                             ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ) ) ) +
                   2 * ( ( 2 * ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) - ( px * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
                       ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                  pz * ( z - zp ) ) *
                                ( z - zp ) ) /
                                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                    pow( z - zp, 2 ) ) ) ) /
                     ( 2. *
                       sqrt( pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) -
                             pow( pz, 2 ) +
                             pow( px - ( ( x - xp ) *
                                         ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                                           ( pow( x - xp, 2 ) +
                                             pow( y - yp, 2 ) +
                                             pow( z - zp, 2 ) ),
                                  2 ) +
                             pow( py - ( ( y - yp ) *
                                         ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                                           ( pow( x - xp, 2 ) +
                                             pow( y - yp, 2 ) +
                                             pow( z - zp, 2 ) ),
                                  2 ) +
                             pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) *
                                         ( z - zp ) ) /
                                           ( pow( x - xp, 2 ) +
                                             pow( y - yp, 2 ) +
                                             pow( z - zp, 2 ) ),
                                  2 ) ) );

  double dMcdy =
      ( 2 *
            ( ( 2 * ( x - xp ) * ( y - yp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) -
              ( py * ( x - xp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ) ) ) *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) +
        2 *
            ( ( 2 * pow( y - yp, 2 ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) -
              ( py * ( y - yp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
              ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                    pow( z - zp, 2 ) ) ) *
            ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) +
        2 *
            ( ( 2 * ( y - yp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) -
              ( py * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ) ) ) *
            ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                     ( z - zp ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) ) /
          ( 2. *
            sqrt(
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) ) +
      ( 2 *
            ( ( 2 * ( x - xp ) * ( y - yp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) -
              ( py * ( x - xp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ) ) ) *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) +
        2 * ( ( 2 * pow( y - yp, 2 ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) - ( py * ( y - yp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) - ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) +
        2 * ( ( 2 * ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) - ( py * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                     ( z - zp ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) ) /
          ( 2. *
            sqrt(
                pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) - pow( pz, 2 ) +
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) );

  double dMcdz =
      ( 2 *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( -( ( pz * ( x - xp ) ) /
                 ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) +
              ( 2 * ( x - xp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 *
            ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( -( ( pz * ( y - yp ) ) /
                 ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) +
              ( 2 * ( y - yp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 *
            ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                     ( z - zp ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( -( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) /
                 ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) -
              ( pz * ( z - zp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) +
              ( 2 * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                pow( z - zp, 2 ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) ) /
          ( 2. *
            sqrt(
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) ) +
      ( 2 *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( -( ( pz * ( x - xp ) ) /
                 ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) +
              ( 2 * ( x - xp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 * ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( -( ( pz * ( y - yp ) ) /
                 ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) +
              ( 2 * ( y - yp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 * ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( -( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) /
                 ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) -
              ( pz * ( z - zp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) +
              ( 2 * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                pow( z - zp, 2 ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) ) /
          ( 2. *
            sqrt(
                pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) - pow( pz, 2 ) +
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) );

  double
      dMcdxp =
          ( 2 *
                ( ( -2 * pow( x - xp, 2 ) *
                    ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                      pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ),
                           2 ) +
                  ( px * ( x - xp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                          pow( z - zp, 2 ) ) +
                  ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) /
                      ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                        pow( z - zp, 2 ) ) ) *
                ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) +
            2 *
                ( ( -2 * ( x - xp ) * ( y - yp ) *
                    ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                      pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ),
                           2 ) +
                  ( px * ( y - yp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                          pow( z - zp, 2 ) ) ) *
                ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) +
            2 *
                ( ( -2 * ( x - xp ) *
                    ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                    ( z - zp ) ) /
                      pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ),
                           2 ) +
                  ( px * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                          pow( z - zp, 2 ) ) ) *
                ( pz -
                  ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                    ( z - zp ) ) /
                      ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                        pow( z - zp, 2 ) ) ) ) /
              ( 2. * sqrt( pow(
                               px - ( ( x - xp ) *
                                      ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                                        ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                          pow( z - zp, 2 ) ),
                               2 ) +
                           pow(
                               py - ( ( y - yp ) *
                                      ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                                        ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                          pow( z - zp, 2 ) ),
                               2 ) +
                           pow(
                               pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) *
                                      ( z - zp ) ) /
                                        ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                          pow( z - zp, 2 ) ),
                               2 ) ) ) +
          ( 2 *
                ( ( -2 * pow( x - xp, 2 ) *
                    ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) /
                      pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                               pow( z - zp, 2 ),
                           2 ) +
                  ( px * ( x - xp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                          pow( z - zp, 2 ) ) +
                  ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) /
                      ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                        pow( z - zp, 2 ) ) ) *
                ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) +
            2 * ( ( -2 * ( x - xp ) * ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) + ( px * ( y - yp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) * ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) + 2 * ( ( -2 * ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) + ( px * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) * ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) ) /
              ( 2. *
                sqrt( pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) -
                      pow( pz, 2 ) +
                      pow( px - ( ( x - xp ) *
                                  ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                                    ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ) ),
                           2 ) +
                      pow( py - ( ( y - yp ) *
                                  ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                                    ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ) ),
                           2 ) +
                      pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) *
                                  ( z - zp ) ) /
                                    ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                      pow( z - zp, 2 ) ),
                           2 ) ) );

  double
      dMcdyp = ( 2 *
                     ( ( -2 * ( x - xp ) * ( y - yp ) *
                         ( px * ( x - xp ) + py * ( y - yp ) +
                           pz * ( z - zp ) ) ) /
                           pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                    pow( z - zp, 2 ),
                                2 ) +
                       ( py * ( x - xp ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) *
                     ( px -
                       ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) +
                 2 *
                     ( ( -2 * pow( y - yp, 2 ) *
                         ( px * ( x - xp ) + py * ( y - yp ) +
                           pz * ( z - zp ) ) ) /
                           pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                    pow( z - zp, 2 ),
                                2 ) +
                       ( py * ( y - yp ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) +
                       ( px * ( x - xp ) + py * ( y - yp ) +
                         pz * ( z - zp ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) *
                     ( py -
                       ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) +
                 2 *
                     ( ( -2 * ( y - yp ) *
                         ( px * ( x - xp ) + py * ( y - yp ) +
                           pz * ( z - zp ) ) *
                         ( z - zp ) ) /
                           pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                    pow( z - zp, 2 ),
                                2 ) +
                       ( py * ( z - zp ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) *
                     ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                pz * ( z - zp ) ) *
                              ( z - zp ) ) /
                                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                  pow( z - zp, 2 ) ) ) ) /
                   ( 2. * sqrt( pow(
                                    px -
                                        ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ),
                                    2 ) +
                                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ), 2 ) + pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ), 2 ) ) ) +
               ( 2 *
                     ( ( -2 * ( x - xp ) * ( y - yp ) *
                         ( px * ( x - xp ) + py * ( y - yp ) +
                           pz * ( z - zp ) ) ) /
                           pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                    pow( z - zp, 2 ),
                                2 ) +
                       ( py * ( x - xp ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) *
                     ( px - ( ( x - xp ) *
                              ( px * ( x - xp ) + py * ( y - yp ) +
                                pz * ( z - zp ) ) ) /
                                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                  pow( z - zp, 2 ) ) ) +
                 2 * ( ( -2 * pow( y - yp, 2 ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) + ( py * ( y - yp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) + ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
                     ( py -
                       ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                        pz * ( z - zp ) ) ) /
                           ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                             pow( z - zp, 2 ) ) ) +
                 2 * ( ( -2 * ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ), 2 ) + ( py * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
                     ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                pz * ( z - zp ) ) *
                              ( z - zp ) ) /
                                ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                  pow( z - zp, 2 ) ) ) ) /
                   ( 2. *
                     sqrt(
                         pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) -
                         pow( pz, 2 ) +
                         pow( px - ( ( x - xp ) *
                                     ( px * ( x - xp ) + py * ( y - yp ) +
                                       pz * ( z - zp ) ) ) /
                                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                         pow( z - zp, 2 ) ),
                              2 ) +
                         pow( py - ( ( y - yp ) *
                                     ( px * ( x - xp ) + py * ( y - yp ) +
                                       pz * ( z - zp ) ) ) /
                                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                         pow( z - zp, 2 ) ),
                              2 ) +
                         pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                                       pz * ( z - zp ) ) *
                                     ( z - zp ) ) /
                                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                         pow( z - zp, 2 ) ),
                              2 ) ) );

  double dMcdzp =
      ( 2 *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( ( pz * ( x - xp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
              ( 2 * ( x - xp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 *
            ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( ( pz * ( y - yp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
              ( 2 * ( y - yp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 *
            ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                     ( z - zp ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) +
              ( pz * ( z - zp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
              ( 2 * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                pow( z - zp, 2 ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) ) /
          ( 2. *
            sqrt(
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) ) +
      ( 2 *
            ( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                    pz * ( z - zp ) ) ) /
                       ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                         pow( z - zp, 2 ) ) ) *
            ( ( pz * ( x - xp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
              ( 2 * ( x - xp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 * ( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( ( pz * ( y - yp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
              ( 2 * ( y - yp ) *
                ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                ( z - zp ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) +
        2 * ( pz - ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) * ( z - zp ) ) / ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) ) *
            ( ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) +
              ( pz * ( z - zp ) ) /
                  ( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ) ) -
              ( 2 * ( px * ( x - xp ) + py * ( y - yp ) + pz * ( z - zp ) ) *
                pow( z - zp, 2 ) ) /
                  pow( pow( x - xp, 2 ) + pow( y - yp, 2 ) + pow( z - zp, 2 ),
                       2 ) ) ) /
          ( 2. *
            sqrt(
                pow( k, 2 ) - pow( px, 2 ) - pow( py, 2 ) - pow( pz, 2 ) +
                pow( px - ( ( x - xp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( py - ( ( y - yp ) * ( px * ( x - xp ) + py * ( y - yp ) +
                                           pz * ( z - zp ) ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) +
                pow( pz - ( ( px * ( x - xp ) + py * ( y - yp ) +
                              pz * ( z - zp ) ) *
                            ( z - zp ) ) /
                              ( pow( x - xp, 2 ) + pow( y - yp, 2 ) +
                                pow( z - zp, 2 ) ),
                     2 ) ) );

  double vertex_errsq =
      cov1( 0, 0 ) * dMcdx * dMcdx + cov1( 1, 1 ) * dMcdy * dMcdy +
      cov1( 2, 2 ) * dMcdz * dMcdz + cov2( 0, 0 ) * dMcdxp * dMcdxp +
      cov2( 1, 1 ) * dMcdyp * dMcdyp + cov2( 2, 2 ) * dMcdzp * dMcdzp +
      2. * dMcdx * dMcdy * cov1( 0, 1 ) + 2. * dMcdx * dMcdz * cov1( 0, 2 ) +
      2. * dMcdy * dMcdz * cov1( 1, 2 ) + 2. * dMcdxp * dMcdyp * cov2( 0, 1 ) +
      2. * dMcdxp * dMcdzp * cov2( 0, 2 ) +
      2. * dMcdyp * dMcdzp * cov2( 1, 2 );

  double momentum_errsq =
      dMcdpx * dMcdpx * cov1( 3, 3 ) + dMcdpy * dMcdpy * cov1( 4, 4 ) +
      dMcdpz * dMcdpz * cov1( 5, 5 ) + dMcdE * dMcdE * cov1( 6, 6 ) +
      // mom v mom cross terms
      +2. * dMcdpx * dMcdpy * cov1( 3, 4 ) +
      2. * dMcdpx * dMcdpz * cov1( 3, 5 ) +
      2. * dMcdpx * dMcdE * cov1( 3, 6 ) +
      2. * dMcdpy * dMcdpz * cov1( 4, 5 ) +
      2. * dMcdpy * dMcdE * cov1( 4, 6 ) + 2. * dMcdpz * dMcdE * cov1( 5, 6 ) +
      // mom vs positon terms
      2 * dMcdx * dMcdpx * cov1( 0, 3 ) + 2. * dMcdy * dMcdpx * cov1( 1, 3 ) +
      2. * dMcdz * dMcdpx * cov1( 2, 3 ) + 2 * dMcdx * dMcdpy * cov1( 0, 4 ) +
      2. * dMcdy * dMcdpy * cov1( 1, 4 ) + 2. * dMcdz * dMcdpy * cov1( 2, 4 ) +
      2 * dMcdx * dMcdpz * cov1( 0, 5 ) + 2. * dMcdy * dMcdpz * cov1( 1, 5 ) +
      2. * dMcdz * dMcdpz * cov1( 2, 5 ) + 2 * dMcdx * dMcdE * cov1( 0, 6 ) +
      2. * dMcdy * dMcdE * cov1( 1, 6 ) + 2. * dMcdz * dMcdE * cov1( 2, 6 );

  std::vector<double> Mcorr_errors;
  Mcorr_errors.push_back( sqrt( vertex_errsq ) );
  Mcorr_errors.push_back( sqrt( vertex_errsq + momentum_errsq ) );

  return Mcorr_errors;
}

double TupleToolSLTools::dpnulowdpperp( double         pperp,
                                        TLorentzVector visible ) {
  double pmag   = visible.P();
  double mBmass = m_Bmass;
  double mvis   = visible.M();
  double mmiss  = sqrt( mBmass * mBmass - mvis * mvis );
  double result =
      ( 2 * pperp +
        ( ( -16 * pmag * pperp * sqrt( 1 - pow( pperp, 2 ) / pow( pmag, 2 ) ) +
            ( -4 * pow( mmiss, 2 ) * pperp + 8 * pow( pperp, 3 ) ) /
                ( pmag * sqrt( 1 - pow( pperp, 2 ) / pow( pmag, 2 ) ) ) +
            ( 16 * pperp *
              ( ( pow( mmiss, 2 ) + pow( mvis, 2 ) ) *
                    ( pow( pmag, 2 ) - 2 * pow( pperp, 2 ) ) +
                pow( mBmass, 2 ) *
                    ( pow( mvis, 2 ) + 2 * pow( pperp, 2 ) ) ) ) /
                sqrt( pow( mmiss, 4 ) * ( pow( mvis, 2 ) + pow( pmag, 2 ) ) +
                      pow( mmiss, 2 ) *
                          ( -4 * pow( pmag, 2 ) * pow( pperp, 2 ) +
                            4 * pow( pperp, 4 ) ) -
                      4 * pow( pperp, 2 ) *
                          ( pow( mvis, 2 ) *
                                ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                            pow( mBmass, 2 ) *
                                ( pow( mvis, 2 ) + pow( pperp, 2 ) ) ) ) ) *
          ( -4 * pmag * ( -pow( mmiss, 2 ) + 2 * pow( pperp, 2 ) ) *
                sqrt( 1 - pow( pperp, 2 ) / pow( pmag, 2 ) ) -
            4 * sqrt( pow( mmiss, 4 ) * ( pow( mvis, 2 ) + pow( pmag, 2 ) ) +
                      pow( mmiss, 2 ) *
                          ( -4 * pow( pmag, 2 ) * pow( pperp, 2 ) +
                            4 * pow( pperp, 4 ) ) -
                      4 * pow( pperp, 2 ) *
                          ( pow( mvis, 2 ) *
                                ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                            pow( mBmass, 2 ) *
                                ( pow( mvis, 2 ) + pow( pperp, 2 ) ) ) ) ) ) /
            ( 32. * pow( pow( mvis, 2 ) + pow( pperp, 2 ), 2 ) ) -
        ( pperp *
          pow(
              -( pow( mmiss, 2 ) * pmag *
                 sqrt( 1 - pow( pperp, 2 ) / pow( pmag, 2 ) ) ) +
                  2 * pmag * pow( pperp, 2 ) *
                      sqrt( 1 - pow( pperp, 2 ) / pow( pmag, 2 ) ) +
                  sqrt( pow( mmiss, 4 ) * ( pow( mvis, 2 ) + pow( pmag, 2 ) ) +
                        pow( mmiss, 2 ) *
                            ( -4 * pow( pmag, 2 ) * pow( pperp, 2 ) +
                              4 * pow( pperp, 4 ) ) -
                        4 * pow( pperp, 2 ) *
                            ( pow( mvis, 2 ) *
                                  ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                              pow( mBmass, 2 ) *
                                  ( pow( mvis, 2 ) + pow( pperp, 2 ) ) ) ),
              2 ) ) /
            pow( pow( mvis, 2 ) + pow( pperp, 2 ), 3 ) ) /
      ( 2. *
        sqrt( pow( pperp, 2 ) +
              pow( -( pow( mmiss, 2 ) * pmag *
                      sqrt( 1 - pow( pperp, 2 ) / pow( pmag, 2 ) ) ) +
                       2 * pmag * pow( pperp, 2 ) *
                           sqrt( 1 - pow( pperp, 2 ) / pow( pmag, 2 ) ) +
                       sqrt( pow( mmiss, 4 ) *
                                 ( pow( mvis, 2 ) + pow( pmag, 2 ) ) +
                             pow( mmiss, 2 ) *
                                 ( -4 * pow( pmag, 2 ) * pow( pperp, 2 ) +
                                   4 * pow( pperp, 4 ) ) -
                             4 * pow( pperp, 2 ) *
                                 ( pow( mvis, 2 ) *
                                       ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                                   pow( mBmass, 2 ) * ( pow( mvis, 2 ) +
                                                        pow( pperp, 2 ) ) ) ),
                   2 ) /
                  ( 4. * pow( pow( mvis, 2 ) + pow( pperp, 2 ), 2 ) ) ) );
  return result;
}

double TupleToolSLTools::dpnuhighdpperp( double         pperp,
                                         TLorentzVector visible ) {
  double pmag   = visible.P();
  double mBmass = m_Bmass;
  double mvis   = visible.M();
  double mmiss  = sqrt( mBmass * mBmass - mvis * mvis );
  double result =
      ( 64 * pperp -
        ( 32 * pperp *
          pow(
              -( pow( mmiss, 2 ) * sqrt( pow( pmag, 2 ) - pow( pperp, 2 ) ) ) +
                  2 * pow( pperp, 2 ) *
                      sqrt( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                  sqrt( pow( pow( mmiss, 2 ) - 2 * pow( pperp, 2 ), 2 ) *
                            ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                        ( pow( mvis, 2 ) + pow( pperp, 2 ) ) *
                            ( -pow( mmiss, 4 ) +
                              4 * pow( pperp, 2 ) *
                                  ( pow( mBmass, 2 ) + pow( pmag, 2 ) -
                                    pow( pperp, 2 ) ) ) ),
              2 ) ) /
            pow( pow( mvis, 2 ) + pow( pperp, 2 ), 3 ) +
        ( ( -4 * sqrt( pow( pmag, 2 ) - pow( pperp, 2 ) ) *
                ( -pow( mmiss, 2 ) + 2 * pow( pperp, 2 ) ) -
            4 * sqrt( pow( pow( mmiss, 2 ) - 2 * pow( pperp, 2 ), 2 ) *
                          ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                      ( pow( mvis, 2 ) + pow( pperp, 2 ) ) *
                          ( -pow( mmiss, 4 ) +
                            4 * pow( pperp, 2 ) *
                                ( pow( mBmass, 2 ) + pow( pmag, 2 ) -
                                  pow( pperp, 2 ) ) ) ) ) *
          ( -16 * pperp * sqrt( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
            ( -4 * pow( mmiss, 2 ) * pperp + 8 * pow( pperp, 3 ) ) /
                sqrt( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
            ( 8 * pperp *
              ( pow( mmiss, 4 ) +
                2 * pow( mmiss, 2 ) *
                    ( pow( pmag, 2 ) - 2 * pow( pperp, 2 ) ) -
                2 * ( 4 * pow( pmag, 2 ) * pow( pperp, 2 ) -
                      6 * pow( pperp, 4 ) +
                      pow( mvis, 2 ) *
                          ( pow( pmag, 2 ) - 2 * pow( pperp, 2 ) ) +
                      pow( mBmass, 2 ) *
                          ( pow( mvis, 2 ) + 2 * pow( pperp, 2 ) ) ) ) ) /
                sqrt( pow( pow( mmiss, 2 ) - 2 * pow( pperp, 2 ), 2 ) *
                          ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                      ( pow( mvis, 2 ) + pow( pperp, 2 ) ) *
                          ( -pow( mmiss, 4 ) +
                            4 * pow( pperp, 2 ) *
                                ( pow( mBmass, 2 ) + pow( pmag, 2 ) -
                                  pow( pperp, 2 ) ) ) ) ) ) /
            pow( pow( mvis, 2 ) + pow( pperp, 2 ), 2 ) ) /
      ( 64. *
        sqrt( pow( pperp, 2 ) +
              pow( -( pow( mmiss, 2 ) *
                      sqrt( pow( pmag, 2 ) - pow( pperp, 2 ) ) ) +
                       2 * pow( pperp, 2 ) *
                           sqrt( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                       sqrt( pow( pow( mmiss, 2 ) - 2 * pow( pperp, 2 ), 2 ) *
                                 ( pow( pmag, 2 ) - pow( pperp, 2 ) ) +
                             ( pow( mvis, 2 ) + pow( pperp, 2 ) ) *
                                 ( -pow( mmiss, 4 ) +
                                   4 * pow( pperp, 2 ) *
                                       ( pow( mBmass, 2 ) + pow( pmag, 2 ) -
                                         pow( pperp, 2 ) ) ) ),
                   2 ) /
                  ( 4. * pow( pow( mvis, 2 ) + pow( pperp, 2 ), 2 ) ) ) );
  return result;
}

double TupleToolSLTools::dq2dpnu( double pnu, TLorentzVector lepton,
                                  TVector3 flight, TLorentzVector visible ) {
  double leptone = lepton.E();
  double leptonx = lepton.Px();
  double leptony = lepton.Py();
  double leptonz = lepton.Pz();
  double flightx = flight.X();
  double flighty = flight.Y();
  double flightz = flight.Z();
  double mvis    = visible.M();
  double px      = visible.Px();
  double py      = visible.Py();
  double pz      = visible.Pz();
  double pmag    = visible.P();
  double mBmass  = m_Bmass;
  double result =
      2 *
      ( leptone + pnu -
        ( flightx * ( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu ) *
          ( leptonx +
            flightx * sqrt( -pow( mBmass, 2 ) +
                            pow( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu,
                                 2 ) ) -
            px ) ) /
            sqrt( -pow( mBmass, 2 ) +
                  pow( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu, 2 ) ) -
        ( flighty * ( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu ) *
          ( leptony +
            flighty * sqrt( -pow( mBmass, 2 ) +
                            pow( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu,
                                 2 ) ) -
            py ) ) /
            sqrt( -pow( mBmass, 2 ) +
                  pow( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu, 2 ) ) -
        ( flightz * ( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu ) *
          ( leptonz +
            flightz * sqrt( -pow( mBmass, 2 ) +
                            pow( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu,
                                 2 ) ) -
            pz ) ) /
            sqrt( -pow( mBmass, 2 ) +
                  pow( sqrt( pow( mvis, 2 ) + pow( pmag, 2 ) ) + pnu, 2 ) ) );

  return result;
}

std::vector<double> TupleToolSLTools::q2_errors( TVector3 v1, TVector3 v2,
                                                 TLorentzVector      p,
                                                 TLorentzVector      lepton,
                                                 Gaudi::SymMatrix7x7 cov1,
                                                 Gaudi::SymMatrix3x3 cov2 ) {
  std::vector<double> q2_err;
  double              x  = v1.Px();
  double              y  = v1.Py();
  double              z  = v1.Pz();
  double              xp = v2.Px();
  double              yp = v2.Py();
  double              zp = v2.Pz();
  double              px = p.Px();
  double              py = p.Py();
  double              pz = p.Pz();
  // double k = p.E();

  TVector3                    flight  = ( v1 - v2 ).Unit();
  std::vector<TLorentzVector> nu_sols = recoNu( p, flight, m_Bmass );
  // std::cout << "number of neutrino solutions " << nu_sols.size() <<
  // std::endl;
  if ( nu_sols.size() < 1 ) {
    return q2_err;
  }
  double pnu_low         = nu_sols[0].P();
  double pnu_high        = nu_sols[1].P();
  double pperp           = ( p.Vect() ).Perp( flight );
  double dpnudpperp_low  = dpnulowdpperp( pperp, p );
  double dpnudpperp_high = dpnuhighdpperp( pperp, p );
  double dq2dnumom_low   = dq2dpnu( pnu_low, lepton, flight, p );
  double dq2dnumom_high  = dq2dpnu( pnu_high, lepton, flight, p );

  // magnificent spaghetti!!!
  double dpperpdx =
      ( -( ( ( 2 * x - 2 * xp ) *
             ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( y - yp, 2 ) ) -
               2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
               pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( z - zp, 2 ) ) +
               pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                pow( z - zp, 2 ) ) -
               2 * px * pz * ( x - xp ) * ( z - zp ) ) ) /
           pow( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                    2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                    pow( zp, 2 ),
                2 ) ) +
        ( pow( py, 2 ) * ( 2 * x - 2 * xp ) +
          pow( pz, 2 ) * ( 2 * x - 2 * xp ) - 2 * px * py * ( y - yp ) -
          2 * px * pz * ( z - zp ) ) /
            ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
              2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
              pow( zp, 2 ) ) ) /
      ( 2. *
        sqrt( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( y - yp, 2 ) ) -
                2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
                pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( z - zp, 2 ) ) +
                pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                 pow( z - zp, 2 ) ) -
                2 * px * pz * ( x - xp ) * ( z - zp ) ) /
              ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                pow( zp, 2 ) ) ) );

  double dpperpdy =
      ( -( ( ( 2 * y - 2 * yp ) *
             ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( y - yp, 2 ) ) -
               2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
               pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( z - zp, 2 ) ) +
               pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                pow( z - zp, 2 ) ) -
               2 * px * pz * ( x - xp ) * ( z - zp ) ) ) /
           pow( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                    2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                    pow( zp, 2 ),
                2 ) ) +
        ( pow( px, 2 ) * ( 2 * y - 2 * yp ) + 2 * pow( pz, 2 ) * ( y - yp ) -
          2 * py * ( px * ( x - xp ) + pz * ( z - zp ) ) ) /
            ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
              2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
              pow( zp, 2 ) ) ) /
      ( 2. *
        sqrt( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( y - yp, 2 ) ) -
                2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
                pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( z - zp, 2 ) ) +
                pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                 pow( z - zp, 2 ) ) -
                2 * px * pz * ( x - xp ) * ( z - zp ) ) /
              ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                pow( zp, 2 ) ) ) );

  double dpperpdz =
      ( -( ( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( y - yp, 2 ) ) -
               2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
               pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( z - zp, 2 ) ) +
               pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                pow( z - zp, 2 ) ) -
               2 * px * pz * ( x - xp ) * ( z - zp ) ) *
             ( 2 * z - 2 * zp ) ) /
           pow( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                    2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                    pow( zp, 2 ),
                2 ) ) +
        ( -2 * px * pz * ( x - xp ) - 2 * py * pz * ( y - yp ) +
          2 * pow( px, 2 ) * ( z - zp ) + 2 * pow( py, 2 ) * ( z - zp ) ) /
            ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
              2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
              pow( zp, 2 ) ) ) /
      ( 2. *
        sqrt( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( y - yp, 2 ) ) -
                2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
                pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( z - zp, 2 ) ) +
                pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                 pow( z - zp, 2 ) ) -
                2 * px * pz * ( x - xp ) * ( z - zp ) ) /
              ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                pow( zp, 2 ) ) ) );
  double dpperpdxp =
      ( -( ( ( -2 * x + 2 * xp ) *
             ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( y - yp, 2 ) ) -
               2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
               pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( z - zp, 2 ) ) +
               pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                pow( z - zp, 2 ) ) -
               2 * px * pz * ( x - xp ) * ( z - zp ) ) ) /
           pow( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                    2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                    pow( zp, 2 ),
                2 ) ) +
        ( pow( py, 2 ) * ( -2 * x + 2 * xp ) +
          pow( pz, 2 ) * ( -2 * x + 2 * xp ) + 2 * px * py * ( y - yp ) +
          2 * px * pz * ( z - zp ) ) /
            ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
              2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
              pow( zp, 2 ) ) ) /
      ( 2. *
        sqrt( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( y - yp, 2 ) ) -
                2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
                pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( z - zp, 2 ) ) +
                pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                 pow( z - zp, 2 ) ) -
                2 * px * pz * ( x - xp ) * ( z - zp ) ) /
              ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                pow( zp, 2 ) ) ) );

  double dpperpdyp =
      ( -( ( ( -2 * y + 2 * yp ) *
             ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( y - yp, 2 ) ) -
               2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
               pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( z - zp, 2 ) ) +
               pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                pow( z - zp, 2 ) ) -
               2 * px * pz * ( x - xp ) * ( z - zp ) ) ) /
           pow( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                    2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                    pow( zp, 2 ),
                2 ) ) +
        ( -2 * pow( pz, 2 ) * ( y - yp ) + pow( px, 2 ) * ( -2 * y + 2 * yp ) +
          2 * py * ( px * ( x - xp ) + pz * ( z - zp ) ) ) /
            ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
              2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
              pow( zp, 2 ) ) ) /
      ( 2. *
        sqrt( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( y - yp, 2 ) ) -
                2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
                pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( z - zp, 2 ) ) +
                pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                 pow( z - zp, 2 ) ) -
                2 * px * pz * ( x - xp ) * ( z - zp ) ) /
              ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                pow( zp, 2 ) ) ) );

  double dpperpdzp =
      ( -( ( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( y - yp, 2 ) ) -
               2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
               pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                pow( z - zp, 2 ) ) +
               pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                pow( z - zp, 2 ) ) -
               2 * px * pz * ( x - xp ) * ( z - zp ) ) *
             ( -2 * z + 2 * zp ) ) /
           pow( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                    2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                    pow( zp, 2 ),
                2 ) ) +
        ( 2 * px * pz * ( x - xp ) + 2 * py * pz * ( y - yp ) -
          2 * pow( px, 2 ) * ( z - zp ) - 2 * pow( py, 2 ) * ( z - zp ) ) /
            ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
              2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
              pow( zp, 2 ) ) ) /
      ( 2. *
        sqrt( ( pow( pz, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( y - yp, 2 ) ) -
                2 * py * ( y - yp ) * ( px * ( x - xp ) + pz * ( z - zp ) ) +
                pow( py, 2 ) * ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) +
                                 pow( z - zp, 2 ) ) +
                pow( px, 2 ) * ( pow( y, 2 ) - 2 * y * yp + pow( yp, 2 ) +
                                 pow( z - zp, 2 ) ) -
                2 * px * pz * ( x - xp ) * ( z - zp ) ) /
              ( pow( x, 2 ) - 2 * x * xp + pow( xp, 2 ) + pow( y, 2 ) -
                2 * y * yp + pow( yp, 2 ) + pow( z, 2 ) - 2 * z * zp +
                pow( zp, 2 ) ) ) );

  // chain rule
  double dq2lowdx   = dq2dnumom_low * dpnudpperp_low * dpperpdx;
  double dq2highdx  = dq2dnumom_high * dpnudpperp_high * dpperpdx;
  double dq2lowdy   = dq2dnumom_low * dpnudpperp_low * dpperpdy;
  double dq2highdy  = dq2dnumom_high * dpnudpperp_high * dpperpdy;
  double dq2lowdz   = dq2dnumom_low * dpnudpperp_low * dpperpdz;
  double dq2highdz  = dq2dnumom_high * dpnudpperp_high * dpperpdz;
  double dq2lowdxp  = dq2dnumom_low * dpnudpperp_low * dpperpdxp;
  double dq2highdxp = dq2dnumom_high * dpnudpperp_high * dpperpdxp;
  double dq2lowdyp  = dq2dnumom_low * dpnudpperp_low * dpperpdyp;
  double dq2highdyp = dq2dnumom_high * dpnudpperp_high * dpperpdyp;
  double dq2lowdzp  = dq2dnumom_low * dpnudpperp_low * dpperpdzp;
  double dq2highdzp = dq2dnumom_high * dpnudpperp_high * dpperpdzp;

  double q2err_low = cov1( 0, 0 ) * dq2lowdx * dq2lowdx +
                     cov1( 1, 1 ) * dq2lowdy * dq2lowdy +
                     cov1( 2, 2 ) * dq2lowdz * dq2lowdz +
                     cov2( 0, 0 ) * dq2lowdxp * dq2lowdxp +
                     cov2( 1, 1 ) * dq2lowdyp * dq2lowdyp +
                     cov2( 2, 2 ) * dq2lowdzp * dq2lowdzp +
                     2. * dq2lowdx * dq2lowdy * cov1( 0, 1 ) +
                     2. * dq2lowdx * dq2lowdz * cov1( 0, 2 ) +
                     2. * dq2lowdy * dq2lowdz * cov1( 1, 2 ) +
                     2. * dq2lowdxp * dq2lowdyp * cov2( 0, 1 ) +
                     2. * dq2lowdxp * dq2lowdzp * cov2( 0, 2 ) +
                     2. * dq2lowdyp * dq2lowdzp * cov2( 1, 2 );

  double q2err_high = cov1( 0, 0 ) * dq2highdx * dq2highdx +
                      cov1( 1, 1 ) * dq2highdy * dq2highdy +
                      cov1( 2, 2 ) * dq2highdz * dq2highdz +
                      cov2( 0, 0 ) * dq2highdxp * dq2highdxp +
                      cov2( 1, 1 ) * dq2highdyp * dq2highdyp +
                      cov2( 2, 2 ) * dq2highdzp * dq2highdzp +
                      2. * dq2highdx * dq2highdy * cov1( 0, 1 ) +
                      2. * dq2highdx * dq2highdz * cov1( 0, 2 ) +
                      2. * dq2highdy * dq2highdz * cov1( 1, 2 ) +
                      2. * dq2highdxp * dq2highdyp * cov2( 0, 1 ) +
                      2. * dq2highdxp * dq2highdzp * cov2( 0, 2 ) +
                      2. * dq2highdyp * dq2highdzp * cov2( 1, 2 );
  q2_err.push_back( sqrt( q2err_low ) );
  q2_err.push_back( sqrt( q2err_high ) );

  return q2_err;
}
