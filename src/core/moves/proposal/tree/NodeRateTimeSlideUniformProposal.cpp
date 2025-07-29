#include <cstddef>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

#include "DistributionUniform.h"
#include "NodeRateTimeSlideUniformProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Proposal.h"
#include "RbVector.h"
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
NodeRateTimeSlideUniformProposal::NodeRateTimeSlideUniformProposal( StochasticNode<Tree> *n, StochasticNode< RbVector<double> > *r ) : Proposal(),
    rates_node( r ),
    variable( n )
{
    // tell the base class to add the node
    addNode( rates_node );
    addNode( variable );
}


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
NodeRateTimeSlideUniformProposal::NodeRateTimeSlideUniformProposal( StochasticNode<Tree> *n, std::vector<StochasticNode<double>* > r ) : Proposal(),
    rates_node( NULL ),
    rates_vector( r ),
    variable( n )
{
    // tell the base class to add the node
    addNode( variable );
    
    for(size_t i = 0; i < rates_vector.size(); i++)
    {
        addNode( rates_vector[i] );
    }

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NodeRateTimeSlideUniformProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NodeRateTimeSlideUniformProposal* NodeRateTimeSlideUniformProposal::clone( void ) const
{
    
    return new NodeRateTimeSlideUniformProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NodeRateTimeSlideUniformProposal::getProposalName( void ) const
{
    static std::string name = "NodeRateTimeSlideUniform";
    
    return name;
}


double NodeRateTimeSlideUniformProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}

/**
 * Perform the proposal.
 *
 * A NodeRateTimeSlide proposal randomly proposes a new age for an internal node
 * uniformly between the age of its parent and its oldest child. The branch
 * rates are then deterministically transformed to preserve constant branch lengths.
 *
 * \return The hastings ratio.
 */
double NodeRateTimeSlideUniformProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    if (tau.getNumberOfTips() <= 2)
    {
        return 0.0;
    }
    
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
        std::cerr << "mvNodeRateTimeSlideUniform has no effect; the tree only contains the root, tips, and sampled ancestors." << std::endl;
        return RbConstants::Double::neginf;
    }

    // pick a random node which is not the root, a tip, or the parent of a SA
    double u = rng->uniform01();
    size_t node_idx = size_t( std::floor(eligible_nodes.size() * u) );
    TopologyNode* node = eligible_nodes[node_idx];
    
    TopologyNode& parent = node->getParent();
    TopologyNode& childA = node->getChild( 0 );
    TopologyNode& childB = node->getChild( 1 );
    
    size_t node_index   = node->getIndex();
    size_t childA_index = childA.getIndex();
    size_t childB_index = childB.getIndex();

    // we need to work with the times
    double parent_age = parent.getAge();
    double my_age     = node->getAge();
    double childA_age = childA.getAge();
    double childB_age = childB.getAge();
    double child_Age  = std::max(childA_age, childB_age);
    
    // now we store all necessary values
    storedNode = node;
    storedAge = my_age;
    
    double prev_node_time   = parent_age - my_age;
    double prev_childA_time = my_age - childA_age;
    double prev_childB_time = my_age - childB_age;

    // draw new ages and compute the hastings ratio at the same time
    double my_new_age = (parent_age-child_Age) * rng->uniform01() + child_Age;
    
    double node_time   = parent_age - my_new_age;
    double childA_time = my_new_age - childA_age;
    double childB_time = my_new_age - childB_age;


    // set the age
    tau.getNode(node_index).setAge( my_new_age );
    
    // set the rates
    if( rates_node != NULL )
    {
        stored_rates = rates_node->getValue();

        std::vector<double>& rates = rates_node->getValue();

        rates[node_index]   = prev_node_time   * rates[node_index]   / node_time;
        rates[childA_index] = prev_childA_time * rates[childA_index] / childA_time;
        rates[childB_index] = prev_childB_time * rates[childB_index] / childB_time;
    }
    else
    {
        stored_rates.clear();

        for(size_t i = 0; i < rates_vector.size(); i++)
        {
            stored_rates.push_back( rates_vector[i]->getValue() );
        }

        rates_vector[node_index]->getValue()   = prev_node_time   * rates_vector[node_index]->getValue()   / node_time;
        rates_vector[childA_index]->getValue() = prev_childA_time * rates_vector[childA_index]->getValue() / childA_time;
        rates_vector[childB_index]->getValue() = prev_childB_time * rates_vector[childB_index]->getValue() / childB_time;
    }

    double ln_numerator   = log(prev_node_time) + log(prev_childA_time) + log(prev_childB_time);
    double ln_denominator = log(node_time)      + log(childA_time)      + log(childB_time);

    return ln_numerator - ln_denominator;
    
}


/**
 *
 */
void NodeRateTimeSlideUniformProposal::prepareProposal( void )
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
void NodeRateTimeSlideUniformProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NodeRateTimeSlideUniformProposal::undoProposal( void )
{
    
    // undo the proposal
    storedNode->setAge( storedAge );
    
    if( rates_node != NULL )
    {
        rates_node->getValue() = stored_rates;
    }
    else
    {
        for(size_t i = 0; i < rates_vector.size(); i++)
        {
            rates_vector[i]->getValue() = stored_rates[i];
        }
    }
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NodeRateTimeSlideUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if( oldN == variable )
    {
        variable = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else if( oldN == rates_node )
    {
        rates_node = static_cast<StochasticNode<RbVector<double> >* >( newN );
    }
    else if( rates_vector.empty() == false )
    {
        for(size_t i = 0; i < rates_vector.size(); i++)
        {
            if(oldN == rates_vector[i])
            {
                rates_vector[i] = static_cast<StochasticNode<double>* >( newN );
            }
        }
    }
    
}


void NodeRateTimeSlideUniformProposal::setProposalTuningParameter(double tp)
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
void NodeRateTimeSlideUniformProposal::tune( double rate )
{
    
}

