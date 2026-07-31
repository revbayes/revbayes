#include "GibbsBuddingTopologyProposal.h"

#include <ostream>

#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RbConstants.h"
#include "RbException.h"
#include "Tree.h"

using namespace RevBayesCore;


GibbsBuddingTopologyProposal::GibbsBuddingTopologyProposal( StochasticNode<Tree> *n ) : Proposal(),
    variable( n ),
    failed( false )
{
    addNode( variable );

    // the uniform draw is the exact conditional only under pure budding
    FossilizedBirthDeathSpeciationProcess* dist = dynamic_cast<FossilizedBirthDeathSpeciationProcess* >( &variable->getDistribution() );

    if ( dist != NULL && dist->hasAnagenesis() == true )
    {
        throw RbException("mvGibbsBuddingTopology is a Gibbs step only under pure budding, so it cannot be used with lambda_a > 0. Use the MH topology moves instead.");
    }

    if ( dist != NULL && dist->hasSymmetricSpeciation() == true )
    {
        throw RbException("mvGibbsBuddingTopology is a Gibbs step only under pure budding, so it cannot be used with beta > 0. Use the MH topology moves instead.");
    }
}


void GibbsBuddingTopologyProposal::cleanProposal( void )
{
}


GibbsBuddingTopologyProposal* GibbsBuddingTopologyProposal::clone( void ) const
{
    return new GibbsBuddingTopologyProposal( *this );
}


const std::string& GibbsBuddingTopologyProposal::getProposalName( void ) const
{
    static std::string name = "GibbsBuddingTopology";

    return name;
}


double GibbsBuddingTopologyProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Draw a whole budding topology from its conditional, which is uniform over the trees the ranges
 * admit. The density is the same for all of them, so the ratio is one and this returns zero.
 */
double GibbsBuddingTopologyProposal::doProposal( void )
{
    FossilizedBirthDeathSpeciationProcess* dist = dynamic_cast<FossilizedBirthDeathSpeciationProcess* >( &variable->getDistribution() );

    if ( dist == NULL )
    {
        throw RbException("mvGibbsBuddingTopology only works on a fossilized birth death speciation process.");
    }

    failed = false;

    stored_tree = variable->getValue();

    if ( dist->redrawTopology() == false )
    {
        failed = true;
    }

    return 0.0;
}


void GibbsBuddingTopologyProposal::prepareProposal( void )
{
}


void GibbsBuddingTopologyProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void GibbsBuddingTopologyProposal::undoProposal( void )
{
    if ( failed == false )
    {
        variable->setValue( stored_tree.clone() );
    }
}


void GibbsBuddingTopologyProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == variable )
    {
        variable = static_cast<StochasticNode<Tree>* >(newN);
    }
}


void GibbsBuddingTopologyProposal::setProposalTuningParameter(double tp)
{
}


void GibbsBuddingTopologyProposal::tune( double rate )
{
}
