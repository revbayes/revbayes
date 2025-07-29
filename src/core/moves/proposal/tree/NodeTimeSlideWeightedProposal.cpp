#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>

#include "DistributionBeta.h"
#include "NodeTimeSlideWeightedProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "DagNode.h"
#include "Proposal.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
NodeTimeSlideWeightedProposal::NodeTimeSlideWeightedProposal( StochasticNode<Tree> *n, size_t b ) : Proposal(),
    variable( n ),
    blocks( b )
{
    // tell the base class to add the node
    addNode( variable );
    
    for (size_t i = 1; i <= blocks; ++i)
    {
        double x = i / (1.0 + blocks);
        double q = RbStatistics::Beta::quantile( 0.25, 0.25, x);
        interval.push_back( q );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NodeTimeSlideWeightedProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NodeTimeSlideWeightedProposal* NodeTimeSlideWeightedProposal::clone( void ) const
{
    
    return new NodeTimeSlideWeightedProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NodeTimeSlideWeightedProposal::getProposalName( void ) const
{
    static std::string name = "NodeTimeSlideWeighted";
    
    return name;
}


double NodeTimeSlideWeightedProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Weighted-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Weighted(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double NodeTimeSlideWeightedProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    // check that there is at least one node which is not the root, a tip, or the parent of a SA
    std::vector<TopologyNode*> eligible_nodes;
    
    for (auto& to_check: tau.getNodes())
    {
        if ( !to_check->isRoot() && !to_check->isTip() && !to_check->isSampledAncestorParent() )
        {
            eligible_nodes.push_back( to_check );
        }
    }
    
    if (eligible_nodes.size() == 0)
    {
        std::cerr << "mvNodeTimeSlide has no effect; the tree only contains the root, tips, and sampled ancestors." << std::endl;
        return RbConstants::Double::neginf;
    }
    
    // pick a random node which is not the root, a tip, or the parent of a SA
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() || node->isSampledAncestorParent() );
    
    TopologyNode& parent = node->getParent();
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child_Age = node->getChild( 1 ).getAge();
    }
    
    // we store all necessary values
    storedNode = node;
    storedAge = my_age;
    
    // approximate the likelihood curve for this node
    std::vector<double> lnl(1,0.0);
    // get the affected dag nodes for the posterior computation
    RbOrderedSet<DagNode*> affected;
    variable->initiateGetAffectedNodes( affected );
    double f = (parent_age - child_Age);
    double marginal = 0.0;
    double prev_x = 0.0;
    double pre_lnl = 0.0;
    for (size_t i = 0; i < blocks; ++i)
    {
        double newAge = interval[i] * f + child_Age;
        tau.getNode(node->getIndex()).setAge( newAge );
        
        double lnLikelihood = variable->getLnProbability();
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            lnLikelihood += (*it)->getLnProbability();
        }
        lnl.push_back( lnLikelihood );
        
        // compute the integral (marginal likelihood)
        marginal += (pre_lnl+lnLikelihood)/2.0 * (interval[i] - prev_x);
        prev_x = interval[i];
        pre_lnl = lnLikelihood;
    }
    // add the final piece of the marginal likelihood
    lnl.push_back( 0.0 );
    marginal += pre_lnl/2.0 * (1.0 - prev_x);
    // normalize the likelihoods
    for (size_t i = 0; i < (blocks+2); ++i) {
        lnl[i] /= marginal;
    }
    
    // randomly draw a new age (using the cdf of the weight function)
    double u = rng->uniform01();
    double proposedAge = 0.0;
    size_t index = 1;
    while ( u > 0 )
    {
        double block = (lnl[index]+lnl[index-1])/2.0 * (interval[index] - interval[index-1]);
        if ( u < block )
        {
            double slope = (lnl[index-1]-lnl[index]) / (interval[index] - interval[index-1]);
            double tmp = sqrt( lnl[index-1]*lnl[index-1] + 2.0*u*slope);
            proposedAge = interval[index-1] + (tmp-lnl[index-1])/slope;
            if ( proposedAge < interval[index-1] || proposedAge > interval[index] )
            {
                proposedAge = interval[index-1] + (-tmp-lnl[index-1])/slope;
            }
        }
        u -= block;
        index++;
    }
    
    // set the age
    tau.getNode(node->getIndex()).setAge( proposedAge );
    
    // compute Hastings ratio (ratio of the weights)
    double weight_old = 1.0, weight_new = 1.0;
    prev_x = 0.0;
    pre_lnl = 0.0;
    bool foundForward = false, foundBackward = false;
    double proposed_x = (proposedAge-child_Age)/f;
    double old_x = (my_age-child_Age)/f;
    for (size_t i = 0; i < blocks; ++i)
    {
        if ( !foundForward && interval[i] > proposed_x)
        {
            foundForward = true;
            weight_new = pre_lnl + (proposed_x-prev_x)/(interval[i]-prev_x)*(lnl[i+1]-pre_lnl);
        }
        if ( !foundBackward && interval[i] > old_x)
        {
            foundBackward = true;
            weight_old = pre_lnl + (old_x-prev_x)/(interval[i]-prev_x)*(lnl[i+1]-pre_lnl);
        }
        if ( foundForward && foundBackward )
        {
            break;
        }
        prev_x = interval[i];
        pre_lnl = lnl[i+1];
    }
    
    return log( weight_old / weight_new );
    
}


/**
 *
 */
void NodeTimeSlideWeightedProposal::prepareProposal( void )
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
void NodeTimeSlideWeightedProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NodeTimeSlideWeightedProposal::undoProposal( void )
{
    
    // undo the proposal
    variable->getValue().getNode(storedNode->getIndex()).setAge( storedAge );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NodeTimeSlideWeightedProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void NodeTimeSlideWeightedProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void NodeTimeSlideWeightedProposal::tune( double rate )
{
    
}

