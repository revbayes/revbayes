#include "ExtendedTipTimeUniformProposal.h"

#include <cstddef>
#include <ostream>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;


ExtendedTipTimeUniformProposal::ExtendedTipTimeUniformProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    node_index( 0 ),
    stored_age( 0.0 ),
    failed( false )
{
    addNode( tree );
}


void ExtendedTipTimeUniformProposal::cleanProposal( void )
{
}


ExtendedTipTimeUniformProposal* ExtendedTipTimeUniformProposal::clone( void ) const
{
    return new ExtendedTipTimeUniformProposal( *this );
}


const std::string& ExtendedTipTimeUniformProposal::getProposalName( void ) const
{
    static std::string name = "ExtendedTipTimeUniform";

    return name;
}


double ExtendedTipTimeUniformProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Draw a new extinction time for a random extinct tip, uniformly on (present, y_i].
 */
double ExtendedTipTimeUniformProposal::doProposal( void )
{
    AbstractFossilizedBirthDeathRangeProcess* dist = dynamic_cast<AbstractFossilizedBirthDeathRangeProcess* >( &tree->getDistribution() );

    if ( dist == NULL )
    {
        throw RbException("mvExtendedTipTimeUniform only works on a fossilized birth death range process.");
    }

    RandomNumberGenerator* rng = GLOBAL_RNG;

    Tree& tau = tree->getValue();

    failed = false;

    // an extant tip is pinned at the present
    std::vector<size_t> tips;
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        const TopologyNode& n = tau.getNode(i);

        if ( n.isTip() == true && n.getTaxon().isExtinct() == true )
        {
            tips.push_back( i );
        }
    }

    if ( tips.empty() == true )
    {
        failed = true;

        return 0.0;
    }

    node_index = tips[ size_t( rng->uniform01() * tips.size() ) ];

    TopologyNode& node = tau.getNode( node_index );

    stored_age = node.getAge();

    // a tip carries its taxon's index, which is how the process reads the tree back
    double present = dist->getPresent();
    double max_age = dist->getMaxExtinctionAge( node_index );

    node.setAge( present + rng->uniform01() * (max_age - present) );

    return 0.0;
}


void ExtendedTipTimeUniformProposal::prepareProposal( void )
{
}


void ExtendedTipTimeUniformProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void ExtendedTipTimeUniformProposal::undoProposal( void )
{
    if ( failed == false )
    {
        tree->getValue().getNode( node_index ).setAge( stored_age );
    }
}


void ExtendedTipTimeUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
}


void ExtendedTipTimeUniformProposal::setProposalTuningParameter(double tp)
{
}


void ExtendedTipTimeUniformProposal::tune( double rate )
{
}
