#include "DistributionBeta.h"
#include "RateAgeProposal.h"
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
RateAgeProposal::RateAgeProposal(StochasticNode<Tree> *n, double a, double p) : Proposal(p),
    tree( n ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( tree );

}


/**
 * Add a variable.
 */
void RateAgeProposal::addRates( std::vector<StochasticNode<double> *> v)
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
void RateAgeProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RateAgeProposal* RateAgeProposal::clone( void ) const
{
    
    return new RateAgeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RateAgeProposal::getProposalName( void ) const
{
    static std::string name = "RateAgeProposal";
    
    return name;
}


double RateAgeProposal::getProposalTuningParameter( void ) const
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
double RateAgeProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    RbOrderedSet<DagNode*> affected;
    tree->getAffectedNodes( affected );
    
    // pick a random node which is not the root, the direct descendant of the root, or the ancestor of a sampled ancestor tip
    TopologyNode* node;
    size_t node_index = 0;
    do {
        double u = rng->uniform01();
        node_index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(node_index);
    } while ( node->isRoot() || node->isTip() || node->isSampledAncestor(true) );
    
    TopologyNode& parent = node->getParent();
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child_Age = node->getChild( 1 ).getAge();
    }
    
    
    // now we store all necessary values
    stored_node = node;
    stored_age  = my_age;
    
    stored_rates.clear();
    for( size_t i = 0; i < branch_rates.size(); ++i )
    {
        
        // make a container for the old values for this vector of rates
        std::vector<double> these_rates;
        
        // add the focal node's rate
        these_rates.push_back( branch_rates[i][node_index]->getValue() );

        // add all of the child rates
        for (size_t j = 0; j < node->getNumberOfChildren(); j++)
        {
            size_t child_index = node->getChild(j).getIndex();
            these_rates.push_back( branch_rates[i][child_index]->getValue() );
        }
        
        // add the vector to the stored rates
        stored_rates.push_back( these_rates );

    }
    
    
    // draw new ages and compute the hastings ratio at the same time
    double m = (my_age-child_Age) / (parent_age-child_Age);
    double a = alpha * m + 1.0;
    double b = alpha * (1.0-m) + 1.0;

    if (a <= 0. || b <= 0. || (!RbMath::isFinite(a) && !RbMath::isFinite(b)))
    {
    	std::cout << a << " -- " << b << " -- " << my_age << " -- " << child_Age << " -- " << parent_age << std::endl;
    }

    double new_m = RbStatistics::Beta::rv(a, b, *rng);
    double my_new_age = (parent_age-child_Age) * new_m + child_Age;
    
    // compute the Hastings ratio
    double forward = RbStatistics::Beta::lnPdf(a, b, new_m);
    double new_a = alpha * new_m + 1.0;
    double new_b = alpha * (1.0-new_m) + 1.0;
    double backward = RbStatistics::Beta::lnPdf(new_a, new_b, m);
    double ln_hastings_ratio = backward - forward;
    
    // set the age
    tau.getNode(node_index).setAge( my_new_age );

    // set the rates
    double pa = node->getParent().getAge();

    
    for (size_t i = 0; i < branch_rates.size(); ++i )
    {
        
        double my_new_rate = (pa - my_age) * stored_rates[i][0] / (pa - my_new_age);
        
        // now we set the new value
        branch_rates[i][node_index]->getValue() = my_new_rate;
        ln_hastings_ratio += log((pa - my_age) / (pa - my_new_age));
        
        for (size_t j = 0; j < node->getNumberOfChildren(); j++)
        {
        
            // compute the new rate for this branch
            size_t child_index = node->getChild(j).getIndex();
            double a = node->getChild(j).getAge();
            double child_new_rate = (my_age - a) * stored_rates[i][j+1] / (my_new_age - a);
            
            // now we set the new value
            branch_rates[i][child_index]->getValue() = child_new_rate;
            
            
            // get the probability ratio of the new rate
            ln_hastings_ratio += log((my_age - a) / (my_new_age - a));
            
        }
        
    }
    
    
    return ln_hastings_ratio;
    
}


/**
 * Prepare the proposal, e.g., pick the UpDown that we want to change.
 * Here we do not need to do any preparation.
 */
void RateAgeProposal::prepareProposal( void )
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
void RateAgeProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void RateAgeProposal::undoProposal( void )
{

    // undo the proposal

    // get this node
    size_t node_index = stored_node->getIndex();

    // reset the age
    tree->getValue().getNode( node_index ).setAge( stored_age );

    // reset the rates
    for( size_t i = 0; i < branch_rates.size(); ++i )
    {
        
        // add the focal node's rate
        branch_rates[i][node_index]->getValue() = stored_rates[i][0];
        
        for ( size_t j = 0; j < stored_node->getNumberOfChildren(); j++ )
        {
            size_t child_index = stored_node->getChild(j).getIndex();
            branch_rates[i][child_index]->getValue() = stored_rates[i][j+1];
        }
        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RateAgeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
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


void RateAgeProposal::setProposalTuningParameter(double tp)
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
void RateAgeProposal::tune( double rate )
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

