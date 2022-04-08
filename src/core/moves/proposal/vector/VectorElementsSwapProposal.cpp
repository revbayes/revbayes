#include "VectorElementsSwapProposal.h"

#include <stddef.h>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
VectorElementsSwapProposal::VectorElementsSwapProposal( std::vector<StochasticNode<double> *> n  ) : Proposal(),
    variables( n )
{
    // tell the base class to add the node
    for (std::vector< StochasticNode<double> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
    {
        addNode( *it );
    }
    
}


void VectorElementsSwapProposal::addIndex( size_t v )
{
    if (v < variables.size())
    {
        indices.insert(v);
    }
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void VectorElementsSwapProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
VectorElementsSwapProposal* VectorElementsSwapProposal::clone( void ) const
{
    return new VectorElementsSwapProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& VectorElementsSwapProposal::getProposalName( void ) const
{
    static std::string name = "VectorElementsSwap";
    
    return name;
}


double VectorElementsSwapProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * lambda)
 * where lambda is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double VectorElementsSwapProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // randomly draw two indices
    size_t chosen_index1;
    size_t chosen_index2;
    if ( indices.size() >= 2 )
    {
        size_t index1 = size_t( floor(rng->uniform01() * double(indices.size())) );
        size_t index2 = size_t( floor(rng->uniform01() * double(indices.size())) );
        while (index1 == index2)
        {
            index2 = size_t( floor(rng->uniform01() * double(indices.size())) );
        }
        
        std::set<size_t>::iterator it1 = indices.begin();
        std::set<size_t>::iterator it2 = indices.begin();
        std::advance(it1, index1);
        std::advance(it2, index2);
        
        chosen_index1 = *it1;
        chosen_index2 = *it2;
    }
    else
    {
        chosen_index1 = size_t( floor(rng->uniform01() * double(variables.size())) );
        chosen_index2 = size_t( floor(rng->uniform01() * double(variables.size())) );
        while (chosen_index1 == chosen_index2)
        {
            chosen_index2 = size_t( floor(rng->uniform01() * double(variables.size())) );
        }
    }
    
    // now we store all necessary values
    stored_chosenidx1 = chosen_index1;
    stored_chosenidx2 = chosen_index2;
    stored_variable1 = variables[chosen_index1]->getValue();
    stored_variable2 = variables[chosen_index2]->getValue();
    
    // swap the values
    variables[chosen_index1]->setValue(new double( stored_variable2 ), false);
    variables[chosen_index2]->setValue(new double( stored_variable1 ), false);

    return 0.0;
    
}


/**
 *
 */
void VectorElementsSwapProposal::prepareProposal( void )
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
void VectorElementsSwapProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void VectorElementsSwapProposal::undoProposal( void )
{
    
    // undo the proposal
    variables[stored_chosenidx1]->setValue(new double( stored_variable1 ), false);
    variables[stored_chosenidx2]->setValue(new double( stored_variable2 ), false);
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void VectorElementsSwapProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
}


void VectorElementsSwapProposal::setProposalTuningParameter(double tp)
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
void VectorElementsSwapProposal::tune( double rate )
{
    
}


