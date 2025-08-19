#include "VectorBinarySwitchProposal.h"


#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
VectorBinarySwitchProposal::VectorBinarySwitchProposal( std::vector< StochasticNode<std::int64_t> *> n, const std::vector<std::int64_t> &i, double sp) : Proposal(),
variables( n ),
indices( i ),
length( indices.size() ),
switch_probability( sp )
{
    // tell the base class to add the node
    // for (std::vector< StochasticNode<std::int64_t> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
    for (size_t i = 0; i < indices.size(); i++)
    {
        size_t index = indices[i]-1;
        addNode( variables[index] );
    }
    
    // if the length is 0 then we use all elements (no indices were provided)
    if ( length == 0 )
    {
        length = variables.size();
    }
    
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void VectorBinarySwitchProposal::cleanProposal( void )
{
    ; //variable->clearTouchedElementIndices();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
VectorBinarySwitchProposal* VectorBinarySwitchProposal::clone( void ) const
{
    
    return new VectorBinarySwitchProposal( *this );
}

/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and Slides the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double VectorBinarySwitchProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // how many elements did we specifically want to update?
    size_t n_indices = indices.size();
    
    // clear update set
    update_set.clear();
    
    // get values to flip
    //    RbVector<std::int64_t> &val = variables->getValue();
    
    // choose one value in the vector to flip
    double u = (size_t)(rng->uniform01() * length);
    size_t first_index = indices[u]-1;
    //    std::int64_t v = variables[u]->getValue() == 0 ? 1 : 0;
    variables[first_index]->setValue( new std::int64_t(variables[first_index]->getValue() == 0 ? 1 : 0) ); //(variables[u]->getValue() == 0 ? 1 : 0) );
    update_set.insert(first_index);
    
    // flip all remaining values in indices according to switch_probability
    for (size_t i = 0; i < length; i++)
    {
        size_t index = ( n_indices == 0 ? i : indices[i]-1 );
        if (rng->uniform01() < switch_probability && index != first_index)
        {
            variables[index]->setValue( new std::int64_t(variables[index]->getValue() == 0 ? 1 : 0) );
            //            val[index] = ( val[index]==0 ? 1 : 0 );
            update_set.insert( index );
        }
    }
    
    return 0.0;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& VectorBinarySwitchProposal::getProposalName( void ) const
{
    static std::string name = "VectorBinarySwitch";
    
    return name;
}


double VectorBinarySwitchProposal::getProposalTuningParameter( void ) const
{
    return switch_probability;
}


std::vector<DagNode*> VectorBinarySwitchProposal::identifyNodesToTouch(void)
{
    std::vector<DagNode*> nodes_to_touch;
 
    for (size_t i = 0; i < indices.size(); i++) {
        size_t index = indices[i]-1;
        if ( update_set.find(index) != update_set.end() ) {
            nodes_to_touch.push_back( nodes[i] );
        }
    }
    return nodes_to_touch;
}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void VectorBinarySwitchProposal::prepareProposal( void )
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
void VectorBinarySwitchProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "switch_probability = ";
    if (name_only == false)
    {
        o << switch_probability;
    }
    
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void VectorBinarySwitchProposal::undoProposal( void )
{
    
    // flip all remaining values in indices according to switch_probability
    for (std::set<size_t>::iterator it = update_set.begin(); it != update_set.end(); it++) {
        size_t index = *it;
        variables[index]->setValue( new std::int64_t(variables[index]->getValue() == 0 ? 1 : 0) );
    }
    
    // clear elements in update_set
    update_set.clear();
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void VectorBinarySwitchProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<std::int64_t> *>(newN);
        }
    }
}


void VectorBinarySwitchProposal::setProposalTuningParameter(double tp)
{
    switch_probability = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void VectorBinarySwitchProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        switch_probability *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        switch_probability /= (2.0 - rate/0.44 );
    }
    
    // bound value s.t. 0 < switch_probability <= 1.0
    if (switch_probability > 1.0)
    {
        switch_probability = 1.0;
    }
    else if (switch_probability < 1.0 / length)
    {
        switch_probability = 1.0 / length;
    }
    
}

