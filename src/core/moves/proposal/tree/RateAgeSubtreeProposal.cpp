#include "DistributionBeta.h"
#include "RateAgeSubtreeProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RateAgeSubtreeProposal::RateAgeSubtreeProposal(StochasticNode<Tree> *n, double a, double p) : Proposal(p),
    tree( n ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( tree );

}


/**
 * Add a variable.
 */
void RateAgeSubtreeProposal::addRates( std::vector<StochasticNode<double> *> v)
{

    // add it to the nodes vector
    for (std::vector<StochasticNode<double>* >::iterator it = v.begin(); it != v.end(); ++it)
    {
        // get the pointer to the current node
        DagNode* the_node = *it;
        addNode( the_node );
    }
    branch_rates.push_back( v );
    
}

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void RateAgeSubtreeProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RateAgeSubtreeProposal* RateAgeSubtreeProposal::clone( void ) const
{
    
    return new RateAgeSubtreeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RateAgeSubtreeProposal::getProposalName( void ) const
{
    static std::string name = "RateAgeSubtreeProposal";
    
    return name;
}


double RateAgeSubtreeProposal::getProposalTuningParameter( void ) const
{
    return alpha;
}


/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( alpha * u )
 * where alpha is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double RateAgeSubtreeProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    RbOrderedSet<DagNode*> affected;
    tree->getAffectedNodes( affected );
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    size_t node_index = 0;
    do {
        double u = rng->uniform01();
        node_index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(node_index);
    } while ( node->isRoot() || node->isTip() );
    
    TopologyNode& parent = node->getParent();
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    
    // now we store all necessary values
    stored_node = node;
    stored_ages = std::vector<double>(tau.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(*node, stored_ages);
    
    // store all the old branch lengths
    std::vector<TopologyNode*> nodes = tau.getNodes();
    size_t num_nodes = nodes.size();
    stored_branch_lengths.clear();
    stored_branch_lengths.resize(num_nodes);
    for(size_t i = 0; i < num_nodes; ++i)
    {
        stored_branch_lengths[i] = tau.getBranchLengthForNode(*nodes[i]);
    }
    
    // lower bound
    double min_age = TreeUtilities::getOldestTipAge(*node);
    
    // draw new ages and compute the hastings ratio at the same time
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

    // reset the rates
    double ln_hastings_ratio = (nNodes > 1 ? log( scaling_factor ) * (nNodes-1) : 0.0 );
    for(size_t j = 0; j < num_nodes; ++j)
    {
        double old_length = stored_branch_lengths[j];
        double new_length = tau.getBranchLengthForNode(*nodes[j]);
        double factor     = old_length / new_length;
        
        if ( new_length != old_length )
        {
            for(size_t i = 0; i < branch_rates.size(); ++i)
            {
                double old_rate = branch_rates[i][j]->getValue();
                double new_rate = old_rate * factor;
                branch_rates[i][j]->getValue() = new_rate;
                ln_hastings_ratio += std::log(factor);
            }
        }
        
    }
    
    return ln_hastings_ratio;
    
}


/**
 * Prepare the proposal, e.g., pick the UpDown that we want to change.
 * Here we do not need to do any preparation.
 */
void RateAgeSubtreeProposal::prepareProposal( void )
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
void RateAgeSubtreeProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void RateAgeSubtreeProposal::undoProposal( void )
{

    // reset the branch rates
    std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
    size_t num_nodes = nodes.size();

    for(size_t j = 0; j < num_nodes; ++j)
    {
        double old_length = stored_branch_lengths[j];
        double new_length = tree->getValue().getBranchLengthForNode(*nodes[j]);
        double factor     = old_length / new_length;
        
        if ( new_length != old_length )
        {
            for(size_t i = 0; i < branch_rates.size(); ++i)
            {
                double old_rate = branch_rates[i][j]->getValue();
                double new_rate = old_rate / factor;
                branch_rates[i][j]->getValue() = new_rate;
            }
        }
        
    }

    // reset the ages
    TreeUtilities::setAges(*stored_node, stored_ages);
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RateAgeSubtreeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else
    {

        for (size_t i = 0; i < branch_rates.size(); ++i)
        {
            
            for (size_t j = 0; j < branch_rates[i].size(); j++)
            {
                if (oldN == branch_rates[i][j])
                {
                    branch_rates[i][j] = static_cast<StochasticNode<double>* >(newN);
                }
            }
            
        }

    }

    
}


void RateAgeSubtreeProposal::setProposalTuningParameter(double tp)
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
void RateAgeSubtreeProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        alpha /= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        alpha *= (2.0 - rate/p);
    }

}

