// local
#include "TupleToolDocas.h"

// from Phys
#include "Kernel/GetIDVAlgorithm.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IDistanceCalculator.h"

using namespace LHCb;
using namespace Gaudi::Units;
using namespace ROOT::Math;

DECLARE_COMPONENT( TupleToolDocas )

/**
 * @brief constructor
 *
 * Properties are three vectors of strings. These should be "aligned": the
 * n-th element of Name names the branch name for the DOCA between the
 * particle the n-th element of Location1 and the n-th element of Location2
 */
TupleToolDocas::TupleToolDocas( const std::string& type,
                                const std::string& name,
                                const IInterface*  parent )
    : TupleToolBase( type, name, parent ),
      m_dva( nullptr ),
      m_dist( nullptr ) {
  declareInterface<IParticleTupleTool>( this );
  declareProperty(
      "Location1", m_locations1,
      "decay descriptor of the particle wrt which the DOCA is computed" );
  declareProperty( "Location2", m_locations2,
                   "decay descriptor of the particle of which the DOCA to the "
                   "reference is computed" );
  declareProperty( "Name", m_name, "name for the branch on the ntuple" );
}

/**
 * @brief clean up memory
 *
 * @return always SUCCESS
 */
StatusCode TupleToolDocas::finalize() {
  // FIXME tried unique pointers, computer said no.
  for ( auto s : m_childSelectors ) {
    delete s.second.first;
    delete s.second.second;
  }
  return TupleToolBase::finalize();
}

/**
 * @brief initialize algorithm
 *
 *  * obtains instance of the distance calculator from the parent
 * DaVinciAlgorithm
 *
 *  * instantiate LoKi::Child::Selectors (which parse the decay descriptor to
 *    obtain descentants from the decay chain
 *
 *  * minimal sanity checks of the properties (name collision and vector
 * misalignment)
 *
 * @return SUCCESS unless there are initialization errors.
 */
StatusCode TupleToolDocas::initialize() {
  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  /// sanity checks for properties
  /// check if vectors have the same length
  if ( m_name.size() != m_locations1.size() ) {
    return Error( "Location1 and Name have different length!",
                  StatusCode::FAILURE );
  }
  if ( m_locations1.size() != m_locations2.size() ) {
    return Error( "Location1 and Location2 have different length!",
                  StatusCode::FAILURE );
  }

  /// quickly test if a branch name is used twice
  for ( auto f = m_name.begin(); f != m_name.end(); ++f ) {
    auto ff = std::find( f + 1, m_name.end(), *f );
    if ( m_name.end() != ff ) {
      return Error( "Two identical branch names in Name",
                    StatusCode::FAILURE );
    }
  }

  /// mostly copy&paste from TupleToolGeometry
  /// in principle we want "the same as LoKi functors", so make sure we use the
  /// LoKi::DistanceCalculator
  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc(), this );
  if ( !m_dva )
    return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE );
  m_dist = m_dva->distanceCalculator( "LoKi::DistanceCalculator" );
  if ( !m_dist )
    return Error( "Unable to retrieve the IDistanceCalculator tool",
                  StatusCode::FAILURE );

  for ( size_t l = 0; l < m_name.size(); ++l ) {
    const auto loc1 = m_locations1[l];
    const auto loc2 = m_locations2[l];
    const auto nam  = m_name[l];
    // FIXME tried unique pointers, computer said no.
    auto selector1 = new LoKi::Child::Selector( loc1 );
    auto selector2 = new LoKi::Child::Selector( loc2 );
    if ( !selector1->valid() ) {
      return Error( std::string( "for " ) + loc1 +
                        " LoKi::Child::Selector returns valid() as false.",
                    StatusCode::FAILURE );
    } else if ( !selector1->valid() ) {
      return Error( std::string( "for " ) + loc2 +
                        " LoKi::Child::Selector returns valid() as false.",
                    StatusCode::FAILURE );
    } else {
      m_childSelectors.emplace( nam, std::make_pair( selector1, selector2 ) );
    }
  }

  return sc;
}

/**
 * @brief compute the docas and fill ntuple
 *
 * fills branch with -999 in case of error
 *
 * @param       unused
 * @param P     the particle for which the DOCA should be evaluated
 * @param head  string addition for branch nameing scheme
 * @param tuple the Tuple from the parent DecayTreeTuple to fill branches
 *
 * @return if tuple filling was successful
 */
StatusCode TupleToolDocas::fill( const Particle*, const Particle* P,
                                 const std::string& head,
                                 Tuples::Tuple&     tuple ) {
  const std::string prefix  = fullName( head ) + "_DOCA_";
  const std::string prefix2 = fullName( head ) + "_DOCACHI2_";
  bool              test    = true;

  // could loop with `for (auto nam: m_name) {` but want to have the index `i`
  // for debug output
  // -> let user know which location failed
  for ( size_t i = 0; i < m_name.size(); ++i ) {
    const auto      nam       = m_name[i];
    const Particle* daughter1 = m_childSelectors[nam].first->child( P );
    const Particle* daughter2 = m_childSelectors[nam].second->child( P );

    // validate that daughters could be found
    if ( nullptr == daughter1 || nullptr == daughter2 ) {
      Error( "didn't find daughters for location", StatusCode::SUCCESS, 20 )
          .ignore();
      if ( nullptr == daughter1 )
        Error( std::string( " failed with location 1: " ) + m_locations1[i],
               StatusCode::SUCCESS, 20 )
            .ignore();
      if ( nullptr == daughter2 )
        Error( std::string( " failed with location 2: " ) + m_locations2[i],
               StatusCode::SUCCESS, 20 )
            .ignore();

      /// fill dummy value -999. due to failure
      test &= tuple->column( prefix + nam, -999. );
      test &= tuple->column( prefix2 + nam, -999. );
    } else {
      double dist = -999.;
      double chi2 = -999.;

      /// actual computation
      StatusCode sc = m_dist->distance( daughter1, daughter2, dist, chi2 );

      /// error handling
      if ( sc.isFailure() ) {
        error() << "doca computation failed" << endmsg;
        error() << "for these two particles" << endmsg;
        error() << " location 1: " << m_locations1[i] << endmsg;
        error() << "particle was " << *daughter1 << endmsg;
        error() << " location 2: " << m_locations2[i] << endmsg;
        error() << "particle was " << *daughter2 << endmsg;

        /// fill dummy value -999. due to failure
        test &= tuple->column( prefix + nam, -999. );
        test &= tuple->column( prefix2 + nam, -999. );
      } else {
        /// computation succeeded, fill result
        test &= tuple->column( prefix + nam, dist );
        test &= tuple->column( prefix2 + nam, chi2 );
      }
    }
    /// inform user upon first filling of branch about branch name
    Info( std::string( "created column " ) + prefix + nam, StatusCode::SUCCESS,
          1 )
        .ignore();
    Info( std::string( "created column " ) + prefix2 + nam,
          StatusCode::SUCCESS, 1 )
        .ignore();
  }

  return StatusCode( test );
}
