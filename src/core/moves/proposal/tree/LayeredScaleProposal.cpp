#include "LayeredScaleProposal.h"

#include <cstddef>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
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
LayeredScaleProposal::LayeredScaleProposal( StochasticNode<Tree> *n, double l, double p ) : Proposal(p),
    variable( n ), lambda( l )
{
    // tell the base class to add the node
    addNode( variable );

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void LayeredScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
LayeredScaleProposal* LayeredScaleProposal::clone( void ) const
{

    return new LayeredScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& LayeredScaleProposal::getProposalName( void ) const
{
    static std::string name = "LayeredScale";

    return name;
}


double LayeredScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 *
 * \return The hastings ratio.
 */
double LayeredScaleProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree& tau = variable->getValue();

    // pick a random node which is not the root or a tip
    TopologyNode* root =  &(tau.getRoot());

    // now we store all necessary values
    storedAges = std::vector<double>(tau.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(*root, storedAges);

    // lower bound
    double min_age = TreeUtilities::getOldestTipAge(*root);

    // upper bound
    double max_age = root->getAge();

    //Below age_layer, all node ages will be rescaled up or down.
    double age_layer = rng->uniform01() * (max_age - min_age) + min_age ;

    std::vector<TopologyNode* > subtrees;
    std::vector<TopologyNode* > nodes = tau.getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (storedAges[i] < age_layer && nodes[i]->getParent().getAge() > age_layer) {
        subtrees.push_back(nodes[i]);
      }
    }

    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );

    size_t nNodes = 0;

    //Now we do the actual rescaling:
    for (size_t i = 0; i < subtrees.size(); ++i) {
      // rescale the subtrees
      TreeUtilities::rescaleSubtree(*subtrees[i], scaling_factor );
      //rescale the branch leading to the subtree
      double parentAge = subtrees[i]->getParent().getAge();
      double currentAge = subtrees[i]->getAge();
      subtrees[i]->setBranchLength(parentAge - age_layer + (age_layer- currentAge)*scaling_factor );
      nNodes+= subtrees[i]->getNumberOfNodesInSubtree(false);
    }

    // compute the Hastings ratio
    double lnHastingsratio = (nNodes > 1 ? log( scaling_factor ) * (nNodes-1) : 0.0 );

    return lnHastingsratio;

}


/**
 *
 */
void LayeredScaleProposal::prepareProposal( void )
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
void LayeredScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void LayeredScaleProposal::undoProposal( void )
{

    // undo the proposal
    TreeUtilities::setAges(variable->getValue().getRoot(), storedAges);

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void LayeredScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    variable = static_cast<StochasticNode<Tree>* >(newN) ;

}


void LayeredScaleProposal::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void LayeredScaleProposal::tune( double rate )
{
  double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        lambda /= (2.0 - rate/p);
    }

}
