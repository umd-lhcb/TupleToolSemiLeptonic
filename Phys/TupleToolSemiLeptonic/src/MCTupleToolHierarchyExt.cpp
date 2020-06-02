// Include files
#include "gsl/gsl_sys.h"
// from Gaudi
#include "GaudiKernel/PhysicalConstants.h"
// local
#include "MCTupleToolHierarchyExt.h"

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/MCParticle.h"
#include "Event/Particle.h"

DECLARE_COMPONENT( MCTupleToolHierarchyExt )

using namespace LHCb;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MCTupleToolHierarchyExt::MCTupleToolHierarchyExt( const std::string& type,
                                                  const std::string& name,
                                                  const IInterface*  parent )
    : TupleToolBase( type, name, parent ) {
  declareInterface<IMCParticleTupleTool>( this );
}

//=============================================================================
StatusCode MCTupleToolHierarchyExt::fill( const LHCb::MCParticle*,
                                          const LHCb::MCParticle* mcp,
                                          const std::string&      head,
                                          Tuples::Tuple&          tuple ) {
  const std::string prefix = fullName( head );
  bool              test   = true;

  int mc_mother_id        = 0;
  int mc_mother_key       = 0;
  int mc_mother_nd        = 0;
  int mc_gd_mother_id     = 0;
  int mc_gd_mother_key    = 0;
  int mc_gd_mother_nd    = 0;
  int mc_gd_gd_mother_id  = 0;
  int mc_gd_gd_mother_key = 0;
  int mc_gd_gd_mother_nd = 0;
  int mc_gd_gd_gd_mother_id = 0;
  int mc_gd_gd_gd_mother_key = 0;
  int mc_gd_gd_gd_mother_nd = 0;

  Gaudi::LorentzVector mc_mother_true_P;
  Gaudi::LorentzVector mc_gd_mother_true_P;
  Gaudi::LorentzVector mc_gd_gd_mother_true_P;
  Gaudi::LorentzVector mc_nu_true_P;

  // pointer is ready, prepare the values:
  if ( mcp ) {
    const MCParticle* mcpmom( 0 );
    mcpmom = mcp->mother();

    if ( mcpmom ) {
      mc_mother_nd = mcpmom->endVertices().front()->products().size();
      mc_mother_id  = mcpmom->particleID().pid();
      mc_mother_key = mcpmom->key();
      mc_mother_true_P = mcpmom->momentum();

      const MCParticle* mcpmom_mom( 0 );
      mcpmom_mom = mcpmom->mother();

      if ( mcp->particleID().pid() == 13 || mcp->particleID().pid() == -13 ) {
        const MCParticle* temp( 0 );
        for ( int counter=0 ; (mcpmom->endVertices().front())->products().size() ; counter++ ) {
          temp = (mcpmom->endVertices()[0])->products()[counter];
          if (temp->particleID().pid() == 14 || temp->particleID().pid() == -14) {
            mc_nu_true_P = temp->momentum();
          }
        }
      }

      if ( mcpmom_mom ) {
        mc_gd_mother_nd = mcpmom_mom->endVertices().front()->products().size();
        mc_gd_mother_id  = mcpmom_mom->particleID().pid();
        mc_gd_mother_key = mcpmom_mom->key();
        mc_gd_mother_true_P = mcpmom_mom->momentum();

        const MCParticle* mcpmom_mom_mom( 0 );
        mcpmom_mom_mom = mcpmom_mom->mother();

        if ( mcpmom_mom_mom ) {
          mc_gd_gd_mother_nd = mcpmom_mom_mom->endVertices().front()->products().size();
          mc_gd_gd_mother_id  = mcpmom_mom_mom->particleID().pid();
          mc_gd_gd_mother_key = mcpmom_mom_mom->key();
          mc_gd_gd_mother_true_P = mcpmom_mom_mom->momentum();

          const MCParticle* mcpmom_mom_mom_mom( 0 );
          mcpmom_mom_mom_mom = mcpmom_mom_mom->mother();

          if ( mcpmom_mom_mom_mom) {
            mc_gd_gd_gd_mother_nd = mcpmom_mom_mom_mom->endVertices().front()->products().size();
            mc_gd_gd_gd_mother_id = mcpmom_mom_mom_mom->particleID().pid();
            mc_gd_gd_gd_mother_key = mcpmom_mom_mom_mom->key();
          }
        }
      }
    }
  }

  // fill the tuple:
  test &= tuple->column( prefix + "_MC_MOTHER_ND", mc_mother_nd );
  test &= tuple->column( prefix + "_MC_MOTHER_ID", mc_mother_id );
  test &= tuple->column( prefix + "_MC_MOTHER_KEY", mc_mother_key );
  test &= tuple->column( prefix + "_MC_MOTHER_TRUEP", mc_mother_true_P );

  test &= tuple->column( prefix + "_MC_GD_MOTHER_ND", mc_gd_mother_nd );
  test &= tuple->column( prefix + "_MC_GD_MOTHER_ID", mc_gd_mother_id );
  test &= tuple->column( prefix + "_MC_GD_MOTHER_KEY", mc_gd_mother_key );
  test &= tuple->column( prefix + "_MC_GD_MOTHER_TRUEP", mc_gd_mother_true_P );

  test &= tuple->column( prefix + "_MC_GD_GD_MOTHER_ND", mc_gd_gd_mother_nd );
  test &= tuple->column( prefix + "_MC_GD_GD_MOTHER_ID", mc_gd_gd_mother_id );
  test &=
      tuple->column( prefix + "_MC_GD_GD_MOTHER_KEY", mc_gd_gd_mother_key );
  test &= tuple->column( prefix + "_MC_GD_GD_MOTHER_TRUEP", mc_gd_gd_mother_true_P );

  test &= tuple->column( prefix + "_MC_GD_GD_GD_MOTHER_ND", mc_gd_gd_gd_mother_nd );
  test &= tuple->column( prefix + "_MC_GD_GD_GD_MOTHER_ID", mc_gd_gd_gd_mother_id );
  test &=
      tuple->column( prefix + "_MC_GD_GD_GD_MOTHER_KEY", mc_gd_gd_gd_mother_key );

  test &= tuple->column( prefix + "_MC_NU_TRUEP", mc_nu_true_P );

  return StatusCode( test );
}
