#include "VectorSimplexSwapProposal.h"

#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "DeterministicNode.h"
//#include "StochasticNode.h"
#include "Cloneable.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
//VectorSimplexSwapProposal::VectorSimplexSwapProposal( StochasticNode< RbVector<Simplex> > *n ) :
VectorSimplexSwapProposal::VectorSimplexSwapProposal( DeterministicNode< RbVector<Simplex> > *n ) :
    Proposal(),
    variable( n ),
    length( variable->getValue().size() ),
    storedIndex1( 0 ),
    storedIndex2( 0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void VectorSimplexSwapProposal::cleanProposal( void )
{
    variable->clearTouchedElementIndices();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
VectorSimplexSwapProposal* VectorSimplexSwapProposal::clone( void ) const
{
    
    return new VectorSimplexSwapProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& VectorSimplexSwapProposal::getProposalName( void ) const
{
    static std::string name = "VectorSimplexSwap";
    
    return name;
}


double VectorSimplexSwapProposal::getProposalTuningParameter( void ) const
{
    return 0.0;
}


/**
 * Perform the proposal.
 *
 * Randomly select 2 indices in the vector and swap the elements.
 *
 * \return The hastings ratio.
 */
double VectorSimplexSwapProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Get vector of simplexes
    RbVector<Simplex> &value = variable->getValue();
    
    // Randomly draw two indices (not necessarily unique)
    this->storedIndex1 = size_t( floor(rng->uniform01() * double(this->length)) );
    this->storedIndex2 = size_t( floor(rng->uniform01() * double(this->length)) );

    // Swap the values located at the chosen indices
    Simplex tmp               = value[this->storedIndex1];
    value[this->storedIndex1] = value[this->storedIndex2];
    value[this->storedIndex2] = tmp;
    
    // This move is symmetric, so lnHastings = 0
    double ln_Hastings_ratio = 0.0;

    return ln_Hastings_ratio;
}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void VectorSimplexSwapProposal::prepareProposal( void )
{
    //-- TODO : Maybe choose random indices here?
}


/**
 * Print the summary of the Proposal.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void VectorSimplexSwapProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    //-- There are no parameters...
    o << " ";
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void VectorSimplexSwapProposal::undoProposal( void )
{
    std::vector<Simplex>& value = variable->getValue();

    // swap the values back
    Simplex tmp               = value[this->storedIndex1];
    value[this->storedIndex1] = value[this->storedIndex2];
    value[this->storedIndex2] = tmp;

    //variable->clearTouchedElementIndices(); //-- TODO : needed?
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void VectorSimplexSwapProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    //variable = static_cast<StochasticNode< RbVector<Simplex> >* >(newN) ;
    variable = static_cast<DeterministicNode< RbVector<Simplex> >* >(newN) ;
    
}


void VectorSimplexSwapProposal::setProposalTuningParameter(double tp)
{
    //-- No tuning parameter for this move
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void VectorSimplexSwapProposal::tune( double rate )
{
    //-- No tuning parameter for this move
}

