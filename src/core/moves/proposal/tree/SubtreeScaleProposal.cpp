#include "DagNodeTypeUtilities.h"
#include "SubtreeScaleProposal.h"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
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

namespace {

// Compute the minimum scale and count exactly the internal ages changed by rescaleSubtree().
// Tips and sampled-ancestor parents are fixed boundaries: the latter's entire subtree is unchanged.
// For an edge from scaled age a to fixed age b, require s*a >= b. Edges with both ends scaled
// preserve their ordering automatically. Take the largest boundary ratio, including ties.
// Descending below a sampled-ancestor parent would incorrectly impose s >= 1 at its fossil tip
// and count unchanged ages in the Jacobian, making an otherwise reversible move one-way.
std::pair<double, size_t> getScalingBoundAndCount(const TopologyNode& node)
{
    double min_scale = 0.0;
    size_t count = 1;
    for (const auto* child : node.getChildren())
    {
        if (child->isTip() || child->isSampledAncestorParent())
        {
            // A boundary at zero imposes no lower bound, even when the parent is also zero.
            if (child->getAge() > 0.0)
                min_scale = std::max(min_scale, child->getAge() / node.getAge());
        }
        else
        {
            auto [child_scale, child_count] = getScalingBoundAndCount(*child);
            min_scale = std::max(min_scale, child_scale);
            count += child_count;
        }
    }
    return {min_scale, count};
}

}

// Draw a new subtree-root age uniformly, and scale the connected internal ages above fixed
// boundaries. The bounds and Hastings term must describe that same set of changed ages.
double SubtreeScaleProposal::doProposal( void )
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
    
    // rescaleSubtree() does nothing at a sampled-ancestor parent. Skip rather than attach a
    // Hastings term to a no-op; a zero root age also cannot define a multiplicative scale.
    if (node->isSampledAncestorParent() || my_age <= 0.0)
        throw RbException(RbException::SKIP_PROPOSAL);

    auto [min_scaling_factor, nNodes] = getScalingBoundAndCount(*node);
    double min_age = min_scaling_factor * my_age;
    if (!(min_age < parent_age))
        throw RbException(RbException::SKIP_PROPOSAL);

    // For selected age r, each boundary contributes r*b/a to this lower bound. The ratio r/a
    // is unchanged by scaling and b is fixed, so the reverse move has the same age interval.
    // The upper bound (the unscaled parent age) and node-selection probability also stay fixed.
    double my_new_age = min_age + (parent_age - min_age) * rng->uniform01();
    double scaling_factor = my_new_age / my_age;

    // Store ages only after all skip checks, and before changing any nodes.
    storedNode = node;
    storedAges = std::vector<double>(tau.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(*node, storedAges);

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

    // The uniform forward/reverse age densities cancel. Of the k scaled ages, the selected
    // age is drawn directly; the other k-1 coordinates are multiplied by s, giving s^(k-1).
    // Count only changed ages: sampled-ancestor parents and everything below them contribute zero.
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

    replaceNodeReference(variable, newN) ;
    
}
