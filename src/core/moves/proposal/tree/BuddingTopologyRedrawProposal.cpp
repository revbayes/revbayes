#include "BuddingTopologyRedrawProposal.h"

#include <ostream>

#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RbConstants.h"
#include "RbException.h"
#include "Tree.h"

using namespace RevBayesCore;


BuddingTopologyRedrawProposal::BuddingTopologyRedrawProposal( StochasticNode<Tree> *n ) : Proposal(),
    variable( n ),
    failed( false )
{
    addNode( variable );

    // the uniform draw is the exact conditional only under pure budding
    FossilizedBirthDeathSpeciationProcess* dist = dynamic_cast<FossilizedBirthDeathSpeciationProcess* >( &variable->getDistribution() );

    if ( dist != NULL && dist->hasAnagenesis() == true )
    {
        throw RbException("mvBuddingTopologyRedraw draws uniformly over compatible topologies, which is the conditional only under pure budding, so lambda_a > 0 is unsupported.");
    }

    if ( dist != NULL && dist->hasSymmetricSpeciation() == true )
    {
        throw RbException("mvBuddingTopologyRedraw draws uniformly over compatible topologies, which is the conditional only under pure budding, so beta > 0 is unsupported.");
    }
}


void BuddingTopologyRedrawProposal::cleanProposal( void )
{
}


BuddingTopologyRedrawProposal* BuddingTopologyRedrawProposal::clone( void ) const
{
    return new BuddingTopologyRedrawProposal( *this );
}


const std::string& BuddingTopologyRedrawProposal::getProposalName( void ) const
{
    static std::string name = "BuddingTopologyRedraw";

    return name;
}


double BuddingTopologyRedrawProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Draw a whole budding topology from its conditional, which is uniform over the trees the ranges
 * admit. The density is the same for all of them, so the ratio is one and this returns zero.
 */
double BuddingTopologyRedrawProposal::doProposal( void )
{
    FossilizedBirthDeathSpeciationProcess* dist = dynamic_cast<FossilizedBirthDeathSpeciationProcess* >( &variable->getDistribution() );

    if ( dist == NULL )
    {
        throw RbException("mvBuddingTopologyRedraw only works on a fossilized birth death speciation process.");
    }

    failed = false;

    stored_tree = variable->getValue();

    if ( dist->redrawTopology() == false )
    {
        failed = true;
    }

    return 0.0;
}


void BuddingTopologyRedrawProposal::prepareProposal( void )
{
}


void BuddingTopologyRedrawProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void BuddingTopologyRedrawProposal::undoProposal( void )
{
    if ( failed == false )
    {
        variable->setValue( stored_tree.clone() );
    }
}


void BuddingTopologyRedrawProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == variable )
    {
        variable = static_cast<StochasticNode<Tree>* >(newN);
    }
}


void BuddingTopologyRedrawProposal::setProposalTuningParameter(double tp)
{
}


void BuddingTopologyRedrawProposal::tune( double rate )
{
}
