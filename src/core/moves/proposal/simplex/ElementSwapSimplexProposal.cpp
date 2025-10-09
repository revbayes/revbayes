#include "ElementSwapSimplexProposal.h"

#include <cstddef>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
ElementSwapSimplexProposal::ElementSwapSimplexProposal( StochasticNode<Simplex> *n ) : SimpleProposal<Simplex>( n ),
    storedValue( Simplex() )
{
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ElementSwapSimplexProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ElementSwapSimplexProposal* ElementSwapSimplexProposal::clone( void ) const
{
    
    return new ElementSwapSimplexProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& ElementSwapSimplexProposal::getProposalName( void ) const
{
    static std::string name = "ElementSwapSimplex";
    
    return name;
}


double ElementSwapSimplexProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * An element swap simplex proposal simply selects two random elements of a simplex
 * and swaps them.
 *
 * \return The hastings ratio.
 */
double ElementSwapSimplexProposal::propose( Simplex &value )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // store the value
    storedValue = value;
    
    // we need to know the number of categories
    size_t cats = value.size();
    
    // randomly draw two indices
    size_t chosen_index_1 = size_t( floor(rng->uniform01()*double(cats)) );
    size_t chosen_index_2 = size_t( floor(rng->uniform01()*double(cats)) );
    while (chosen_index_1 == chosen_index_2)
    {
        chosen_index_2 = size_t( floor(rng->uniform01()*double(cats)) );
    }
   
    // swap the values
    double value_1 = value[chosen_index_1];
    double value_2 = value[chosen_index_2];
    value[chosen_index_1] = value_2;
    value[chosen_index_2] = value_1;
    
    
    double ln_Hastins_ratio = 0;
    
    return ln_Hastins_ratio;
}


/**
 *
 */
void ElementSwapSimplexProposal::prepareProposal( void )
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
void ElementSwapSimplexProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void ElementSwapSimplexProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->setValue( storedValue.clone() );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void ElementSwapSimplexProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    variable = static_cast<StochasticNode<Simplex >* >(newN) ;
}


void ElementSwapSimplexProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 */
void ElementSwapSimplexProposal::tune( double rate )
{
    
}

