// $Id: TupleToolMCDaughters.cpp,v 1.17 2010-01-26 15:39:26 rlambert Exp $
// Include files
#include "gsl/gsl_sys.h"
// from Gaudi
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/PhysicalConstants.h"
// local
#include "TupleToolMCDaughters.h"

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/Particle.h"
#include "Event/MCParticle.h"

// kernel
#include "Kernel/IParticle2MCAssociator.h"
using namespace LHCb;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolMCDaughters::TupleToolMCDaughters( const std::string& type,
                                    const std::string& name,
                                    const IInterface* parent )
  : TupleToolBase ( type, name , parent )
  , m_toolList(1,"MCTupleToolKinematic")
{
  declareInterface<IParticleTupleTool>(this);

  // The names of MCTupleTools to use on the associated mcp
  declareProperty( "ToolList", m_toolList  );
  declareProperty( "Mother", m_Mother = false);
  declareProperty( "NMin", m_nMin = 5);
  // MC associators to try, in order
  m_p2mcAssocTypes.push_back( "DaVinciSmartAssociator" );
  m_p2mcAssocTypes.push_back( "MCMatchObjP2MCRelator"  );
  declareProperty( "IP2MCPAssociatorTypes", m_p2mcAssocTypes );
}

//=============================================================================

StatusCode TupleToolMCDaughters::initialize()
{
  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  // the MC associators
  for ( std::vector<std::string>::const_iterator iMCAss = m_p2mcAssocTypes.begin();
        iMCAss != m_p2mcAssocTypes.end(); ++iMCAss )
  {
    m_p2mcAssocs.push_back( tool<IParticle2MCAssociator>(*iMCAss,this) );
  }

  // initialise the tuple tools
  std::sort( m_toolList.begin(), m_toolList.end() );
  m_toolList.erase( std::unique(m_toolList.begin(),m_toolList.end()), m_toolList.end() );
  for ( std::vector<std::string>::const_iterator it = m_toolList.begin();
        m_toolList.end()!=it ; ++it )
  {
    if (msgLevel(MSG::VERBOSE)) verbose() << "Adding the tool " << *it << endmsg ;
    IMCParticleTupleTool* aTool = tool<IMCParticleTupleTool>(*it,this);
    if ( aTool )
    {
      m_mcTools.push_back(aTool);
    }
    else
    {
      Warning("There was a problem retrieving " + *it +" , this tool will be ignored").ignore();
      std::cout << "\n\n\nThere was a problem retrieving " + *it +" , this tool will be ignored\n\n\n" << std::endl;
    }
  }

  if (msgLevel(MSG::VERBOSE))
  {
    verbose() << "Completed TupleTool intialisation, "
              << m_mcTools.size()
              << " tools added " << endmsg ;
  }

  return sc;
}

//=============================================================================

StatusCode TupleToolMCDaughters::fill( const LHCb::Particle*
                                   , const LHCb::Particle* P
                                   , const std::string& head
                                   , Tuples::Tuple& tuple )
{
  const std::string prefix = fullName(head);

  //warning() << "filling" <<endmsg;
  Assert( !m_p2mcAssocs.empty(),
          "The DaVinci smart associator(s) have not been initialized!");

  // int mcPid = 0;
  bool test = true;
  const LHCb::MCParticle* mcpP(NULL);
  const LHCb::MCParticle* mcp(NULL);
  std::vector<Gaudi::LorentzVector> relatedVectors;
  std::vector<int> relatedIDs;
  std::string related = "DAUGHTER"; // if mother not set, the particles retrived are daughters of P; else sisters

  if ( P )
  {
    //assignedPid = P->particleID().pid();
    if (msgLevel(MSG::VERBOSE)) verbose() << "Getting related MCP to " << P << endmsg ;
	  bool foundMCP=false;
    for ( std::vector<IParticle2MCAssociator*>::const_iterator iMCAss = m_p2mcAssocs.begin();
          iMCAss != m_p2mcAssocs.end(); ++iMCAss )
    {
      mcpP = (*iMCAss)->relatedMCP(P);
      if ( mcpP ) {foundMCP=true;break;}
    }
    if(!foundMCP){
      warning() << "NOT FOUND ASSOCIATED MCPARTICLE" << endmsg;
      return StatusCode(true);
    }
    if (msgLevel(MSG::VERBOSE))
      verbose() << "Got mcp " << mcpP << endmsg ;
  }
   
  if(m_Mother) { // if requested, make mcp the mother particle (mcpP not used below)
    mcp = mcpP->mother();
    related = "SISTER";
  } else {
    mcp=mcpP;
  }
  if(m_Mother && !(mcp)) {
    warning() << "NOT FOUND MOTHER " << endmsg;
    return StatusCode(test);
  }
  //warning() << "Mother ID " << mcp->particleID().abspid() << endmsg;
  // pointer is ready, prepare the values:
  const SmartRefVector < LHCb::MCVertex >& EndVertices = mcp->endVertices();
  //if(EndVertices.size() < 1){warning() << "ENDVERTICES SIZE: " << EndVertices.size() << endmsg;return StatusCode(test);}
  //if(EndVertices.size() > 1){warning() << "ENDVERTICES SIZE: " << EndVertices.size() << endmsg;return StatusCode(test);}
  SmartRefVector< MCVertex >::const_iterator itV = EndVertices.begin();
  if(EndVertices.size() > 1) ++ itV;
  const LHCb::MCVertex *endV = (*itV);
  //warning() << "ENDVERTICES SIZE: " << EndVertices.size() << endmsg;
  const SmartRefVector< LHCb::MCParticle > daughters = endV->products();

  int i=0;
  std::stringstream sout;
  // Greg's old code
  // warning() << "looping through" <<endmsg;
  // SmartRefVector< LHCb::MCParticle >::const_iterator itD;
  // for( itD = daughters.begin(); daughters.end()!=itD; ++itD ){
  //   if((*itD)->particleID().abspid() == 22) continue;
  //   //warning() << "Daughter ID " << endmsg;
  //   //warning() << (*itD)->particleID().abspid() << endmsg;
  //   sout << i;
  //   double tid = (double) (*itD)->particleID().abspid();
  //   //test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_ID",  (*itD)->particleID().abspid());
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_ID", tid);
  //   std::cout << tid << std::endl;
  //   //test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_ID", (*itD)->particleID().abspid());
  //   //warning() << "After ID" << endmsg;
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PE", (double) (*itD)->momentum().E());
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PX", (double) (*itD)->momentum().Px());
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PY", (double) (*itD)->momentum().Py());
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PZ", (double) (*itD)->momentum().Pz());
  //   //warning() << "After ALL" << endmsg;
  //   sout.str("");
  //   i++;
  //   //if(i==6) warning() << "REACHED 6" << endmsg;
  //   //if(i>6) warning() << "ABOVE 6" << endmsg;
  //   if(i>m_nMin) {warning() << "REACHED MAXIMUM nDAUGHTERS" << endmsg; break;}
  // }

  // //warning() << "empty filling" <<endmsg;
	// //warning() << "After LOOP" << endmsg;
  // for(;i<m_nMin;i++){
  //   sout << i;
  //   //if(i==2) continue;
  //   //warning() << "HERE?? " << i << endmsg;
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_ID", (double) 0.);
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PE", (double) 0.);
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PX", (double) 0.);
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PY", (double) 0.);
  //   test &= tuple->column( prefix+"_DAUGHTER"+sout.str()+"_PZ", (double) 0.);
  //   sout.str("");
  // }
  // if(m_Mother) {
  //   test &= tuple->column( prefix+"_TESTMOTHER_ID", (double) mcp->particleID().abspid());
  //   test &= tuple->column( prefix+"_TESTMOTHER_PE", (double) mcp->momentum().E());
  //   test &= tuple->column( prefix+"_TESTMOTHER_PX", (double) mcp->momentum().Px());
  //   test &= tuple->column( prefix+"_TESTMOTHER_PY", (double) mcp->momentum().Py());
  //   test &= tuple->column( prefix+"_TESTMOTHER_PZ", (double) mcp->momentum().Pz());
  // }
  // //warning() << "done" <<endmsg;
  SmartRefVector< LHCb::MCParticle >::const_iterator itD;
  for( itD = daughters.begin(); daughters.end()!=itD; ++itD ){
    i++;
    if(i>m_nMin) {warning() << "REACHED MAXIMUM nDAUGHTERS" << endmsg; break;}
    if((*itD)->particleID().abspid() == 22) continue;
    relatedIDs.push_back( (*itD)->particleID().pid() );
    relatedVectors.push_back( (*itD)->momentum() );
  }
  i=0;
  for (;i<m_nMin;i++) {
    sout << i;
    if (i >= relatedIDs.size()) {
      test &= tuple->column( prefix+"_"+related+sout.str()+"_ID", 0);
      test &= tuple->column( prefix+"_"+related+sout.str()+"_P", Gaudi::LorentzVector(0, 0, 0, 0));
    } else {
      test &= tuple->column( prefix+"_"+related+sout.str()+"_ID", relatedIDs[i]);
      test &= tuple->column( prefix+"_"+related+sout.str()+"_P", relatedVectors[i]);
    }
    sout.str("");
    if(m_Mother) { // already know from above that if m_Mother is true, mcp must not be null ptr
      test &= tuple->column( prefix+"_MOM_ID", mcp->particleID().pid());
      test &= tuple->column( prefix+"_MOM_P", mcp->momentum());
    }
  }
  return StatusCode(test);
}

//=============================================================================

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_COMPONENT( TupleToolMCDaughters )
