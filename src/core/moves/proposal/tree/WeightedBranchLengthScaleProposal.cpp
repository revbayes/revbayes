#include "WeightedBranchLengthScaleProposal.h"

#include <cmath>
#include <iostream>

#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
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
WeightedBranchLengthScaleProposal::WeightedBranchLengthScaleProposal( StochasticNode<Tree> *t, size_t n, double a ) : Proposal(),
    tree( t ),
    num_breaks( n ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( tree );
    
    interval = std::vector<double>( num_breaks, 1.0 );
    for (size_t i = 1; i <= num_breaks; ++i)
    {
        double x = i / (1.0 + num_breaks);
        double q = RbStatistics::Beta::quantile( 0.25, 0.25, x);
        interval[i-1] = q;
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void WeightedBranchLengthScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
WeightedBranchLengthScaleProposal* WeightedBranchLengthScaleProposal::clone( void ) const
{
    
    return new WeightedBranchLengthScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& WeightedBranchLengthScaleProposal::getProposalName( void ) const
{
    static std::string name = "Weighted-BranchLengthScale";
    
    return name;
}


double WeightedBranchLengthScaleProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 *
 *
 * \return The hastings ratio.
 */
double WeightedBranchLengthScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree& tau = tree->getValue();


    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node = NULL;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() == true || node->getParent().isRoot() == true );

    // get the parent node
    TopologyNode& parent = node->getParent();
    
    // get the sibling node
    TopologyNode* sibling = &parent.getChild( 0 );
    // check if we got the correct child
    if ( node == sibling )
    {
        sibling = &parent.getChild( 1 );
    }
    
    // now we need the branch length of the parent (towards the grant parent) and of the sibling
    double branch_length_parent  = parent.getBranchLength();
    double branch_length_sibling = sibling->getBranchLength();
    double branch_length_total   = branch_length_parent + branch_length_sibling;

    // store the values
    stored_branch_index = node->getIndex();
    stored_value = branch_length_parent / branch_length_total;
    
    std::vector<double> lnl = std::vector<double>(num_breaks+2, 1.0);
    
    // get the affected dag nodes for the posterior computation
    RbOrderedSet<DagNode*> affected;
    tree->initiateGetAffectedNodes( affected );
    
    for (size_t i = 0; i < num_breaks; ++i)
    {
        // compute the new parent branch length and set it
        double new_parent_bl = interval[i] * branch_length_total;
        parent.setBranchLength( new_parent_bl );
        
        // compute the new sibling branch length and set it
        double new_sibling_bl = (1.0-interval[i]) * branch_length_total;
        sibling->setBranchLength( new_sibling_bl );
        
        // flag the tree as dirty
        tree->touch();
        
        double ln_likelihood = tree->getLnProbability();
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            ln_likelihood += (*it)->getLnProbability();
        }
        lnl[i+1] = ln_likelihood;
    }
    // project the first log likelihood
    lnl[0] = lnl[1] - (lnl[2]-lnl[1])*interval[0]/(interval[1]-interval[0]);
    // project the last log likelihood
    lnl[num_breaks+1] = lnl[num_breaks] + (lnl[num_breaks]-lnl[num_breaks-1])*(1.0-interval[num_breaks-1])/(interval[num_breaks-1]-interval[num_breaks-2]);
    
    // find the maximum lnl
    double max_lnl = lnl[0];
    for ( size_t i=1; i<(num_breaks+2); ++i )
    {
        if ( max_lnl < lnl[i] )
        {
            max_lnl = lnl[i];
        }
    }
    
    // now transform into likelihood
    std::vector<double> likelihoods = std::vector<double>(num_breaks+2, 1.0);
    for ( size_t i=0; i<(num_breaks+2); ++i )
    {
        likelihoods[i] = exp(lnl[i]-max_lnl);
    }

    
    // compute the integral (marginal likelihood)
    double prev_x = 0.0;
    double pre_like = likelihoods[0];
    double marginal = 0;
    for (size_t i = 0; i < num_breaks; ++i)
    {
        marginal += (pre_like+likelihoods[i+1])/2.0 * (interval[i] - prev_x);
        prev_x = interval[i];
        pre_like = likelihoods[i+1];
    }
    marginal += (pre_like+likelihoods[num_breaks+1])/2.0 * (1.0 - prev_x);
    
    // normalize the likelihoods
    for (size_t i = 0; i < (num_breaks+2); ++i)
    {
        likelihoods[i] /= marginal;
    }
    
    // randomly draw a new branch fraction (using the cdf of the weight function)
    double u = rng->uniform01();
    double proposed_branch_fraction = 0.0;
    size_t index = 1;
    while ( u > 0 )
    {
        double x1 = index > 1 ? interval[index-2] : 0.0;
        double x2 = index < (num_breaks+1) ? interval[index-1] : 1.0;
        double y1 = likelihoods[index-1];
        double y2 = likelihoods[index];
        
        // compute the total probability for this interval
        double block = (y1+y2)/2.0 * (x2-x1);
        if ( u < block )
        {
            // compute the slope
            double alpha = (y2-y1) / (x2-x1);
            
            double sol = 0.0;
            if ( alpha == 0.0 )
            {
                sol = u/y1;
            }
            else
            {
                // solve the quadratic equation
                double p = 2.0*y1/alpha;
                double q = -2.0*u/alpha;
                double tmp = sqrt( p*p/4.0 - q );
                double sol1 = -p/2.0 + tmp;
                double sol2 = -p/2.0 - tmp;
                if ( sol1 < ( x2-x1 ) && sol1 > 0 )
                {
                    sol = sol1;
                }
                else
                {
                    sol = sol2;
                }
            }
            
            //
            proposed_branch_fraction = x1 + sol;
            if ( proposed_branch_fraction < x1 || proposed_branch_fraction > x2 )
            {
                throw RbException("Wrong proposal");
            }
        }
        u -= block;
        index++;
    }
    proposed_branch_fraction *= branch_length_total;
    
    // set the branch lengths
    parent.setBranchLength( proposed_branch_fraction );
    sibling->setBranchLength( branch_length_total-proposed_branch_fraction );
    
    // compute Hastings ratio (ratio of the weights)
    double weight_old = 1.0, weight_new = 1.0;
    prev_x = 0.0;
    pre_like = likelihoods[0];
    bool found_forward = false, found_backward = false;
    double proposed_x = proposed_branch_fraction / branch_length_total;
    double old_x = branch_length_parent / branch_length_total;
    for (size_t i = 0; i < num_breaks; ++i)
    {
        if ( !found_forward && interval[i] > proposed_x)
        {
            found_forward = true;
            double alpha = (likelihoods[i+1]-pre_like) / (interval[i]-prev_x);
            weight_new = pre_like + (proposed_x-prev_x)*alpha;
        }
        if ( !found_backward && interval[i] > old_x)
        {
            found_backward = true;
            double alpha = (likelihoods[i+1]-pre_like) / (interval[i]-prev_x);
            weight_old = pre_like + (old_x-prev_x)*alpha;
        }
        if ( found_forward && found_backward )
        {
            break;
        }
        prev_x = interval[i];
        pre_like = likelihoods[i+1];
    }
    
    // if we haven't found it until the last interval
    if ( !found_forward )
    {
        double alpha = (likelihoods[num_breaks+1]-pre_like) / (1.0-prev_x);
        weight_new = pre_like + (proposed_x-prev_x)*alpha;
    }
    if ( !found_backward )
    {
        double alpha = (likelihoods[num_breaks+1]-pre_like) / (1.0-prev_x);
        weight_old = pre_like + (old_x-prev_x)*alpha;
    }
    
    return log( weight_old / weight_new );
}


/**
 *
 */
void WeightedBranchLengthScaleProposal::prepareProposal( void )
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
void WeightedBranchLengthScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void WeightedBranchLengthScaleProposal::undoProposal( void )
{
    
    Tree& tau = tree->getValue();
    
    TopologyNode* node = &tau.getNode(stored_branch_index);
    
    // get the parent node
    TopologyNode& parent = node->getParent();
    
    // get the sibling node
    TopologyNode* sibling = &parent.getChild( 0 );
    // check if we got the correct child
    if ( node == sibling )
    {
        sibling = &parent.getChild( 1 );
    }
    
    // now we need the branch length of the parent (towards the grant parent) and of the sibling
    double branch_length_parent  = parent.getBranchLength();
    double branch_length_sibling = sibling->getBranchLength();
    double branch_length_total   = branch_length_parent + branch_length_sibling;
        
    // compute the new parent branch length and set it
    double old_parent_bl = stored_value * branch_length_total;
    parent.setBranchLength( old_parent_bl, false );
        
    // compute the new sibling branch length and set it
    double old_sibling_bl = (1.0-stored_value) * branch_length_total;
    sibling->setBranchLength( old_sibling_bl, false );
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void WeightedBranchLengthScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
    
}


void WeightedBranchLengthScaleProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void WeightedBranchLengthScaleProposal::tune( double rate )
{
    
//    if ( rate > 0.234 )
//    {
//        alpha *= (1.0 + ((rate-0.234)/0.766) );
//    }
//    else
//    {
//        alpha /= (2.0 - rate/0.234 );
//    }
    
}


