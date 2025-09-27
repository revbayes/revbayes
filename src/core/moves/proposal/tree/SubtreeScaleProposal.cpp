#include "SubtreeScaleProposal.h"

#include <cstddef>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TreeUtilities.h"
#include "Cloneable.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SubtreeScaleProposal::SubtreeScaleProposal( StochasticNode<Tree> *n ) : Proposal(),
    variable( n )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SubtreeScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SubtreeScaleProposal* SubtreeScaleProposal::clone( void ) const
{
    
    return new SubtreeScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SubtreeScaleProposal::getProposalName( void ) const
{
    static std::string name = "SubtreeScale";
    
    return name;
}


double SubtreeScaleProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}

std::pair<double,double> getOldestTipAgeAndScalingFactor(TopologyNode& n)
{
    double min_scaling_factor = 0;
    double max_tip_age = n.getAge();

    if ( not n.isTip() )
    {

        // assertion that we have binary trees
#ifdef ASSERTIONS_TREE
        if ( n->getNumberOfChildren() != 2 )
        {
            throw RbException("Oldest tip is only implemented for binary trees!");
        }
#endif

        auto [left_scale,  max_left_tip_age ] = getOldestTipAgeAndScalingFactor( n.getChild(0) );
        auto [right_scale, max_right_tip_age] = getOldestTipAgeAndScalingFactor( n.getChild(1) );

        max_tip_age = std::max(max_left_tip_age, max_right_tip_age);

        min_scaling_factor = max_tip_age / n.getAge();
        min_scaling_factor = std::max(min_scaling_factor, left_scale);
        min_scaling_factor = std::max(min_scaling_factor, right_scale);
    }

    return {min_scaling_factor, max_tip_age};
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
LogDensity SubtreeScaleProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        auto& node = tau.getNode(i);

        assert(node.isRoot() or node.getAge() <= node.getParent().getAge());
    }

    // pick a random node which is not the root or a tip
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() );
    
    TopologyNode& parent = node->getParent();
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    
    // now we store all necessary values
    storedNode = node;
    storedAges = std::vector<double>(tau.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(*node, storedAges);

    // lower bound
    auto [min_scaling_factor, max_tip_age] = getOldestTipAgeAndScalingFactor(*node);
    
    // draw new ages and compute the hastings ratio at the same time
    double min_age = max_tip_age;
    min_age = std::max(max_tip_age, min_scaling_factor * node->getAge());

    double my_new_age = min_age + (parent_age - min_age) * rng->uniform01();
    
    double scaling_factor = my_new_age / my_age;
    
    size_t nNodes = node->getNumberOfNodesInSubtree(false);
    
    // rescale the subtrees
    TreeUtilities::rescaleSubtree(*node, scaling_factor );
    
    if (min_age != 0.0)
    {
        for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
        {
            if (tau.getNode(i).getAge() < 0.0) {
                return RbConstants::Double::neginf;
            }
        }
    }
    
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        auto& node = tau.getNode(i);

        assert(node.isRoot() or node.getAge() <= node.getParent().getAge());
    }

    // compute the Hastings ratio
    double lnHastingsratio = (nNodes > 1 ? log( scaling_factor ) * (nNodes-1) : 0.0 );
    
    return lnHastingsratio;
    
}


/**
 *
 */
void SubtreeScaleProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void SubtreeScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SubtreeScaleProposal::undoProposal( void )
{
    
    // undo the proposal
    TreeUtilities::setAges(*storedNode, storedAges);

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SubtreeScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}
