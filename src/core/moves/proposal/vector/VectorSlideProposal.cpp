#include "VectorSlideProposal.h"

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
VectorSlideProposal::VectorSlideProposal( StochasticNode< RbVector<double> > *n, const std::vector<std::int64_t> &i, double l) : Proposal(),
    variable( n ),
    indices( i ),
    lambda( l ),
    length( indices.size() ),
    storedSlidingFactor( 0.0 )
{
    // tell the base class to add the node
    addNode( variable );
    
    // if the length is 0 then we use all elements (no indices were provided)
    if ( length == 0 )
    {
        length = variable->getValue().size();
    }
    
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void VectorSlideProposal::cleanProposal( void )
{
    variable->clearTouchedElementIndices();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
VectorSlideProposal* VectorSlideProposal::clone( void ) const
{
    
    return new VectorSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& VectorSlideProposal::getProposalName( void ) const
{
    static std::string name = "VectorSlide";
    
    return name;
}


double VectorSlideProposal::getProposalTuningParameter( void ) const
{
    return lambda;
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
double VectorSlideProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    RbVector<double> &val = variable->getValue();
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double delta  = ( lambda * ( u - 0.5 ) );
    
    // copy value
    storedSlidingFactor = delta;
    
    // how many elements did we specifically want to update?
    size_t n_indices = indices.size();
    
    for (size_t i=0; i<length; ++i)
    {
        size_t index = ( n_indices == 0 ? i : indices[i]);
        val[index] += delta;
    }
    
    return 0.0;
}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void VectorSlideProposal::prepareProposal( void )
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
void VectorSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void VectorSlideProposal::undoProposal( void )
{
    RbVector<double>& v = variable->getValue();
    
    // how many elements did we specifically want to update?
    size_t n_indices = indices.size();
    
    for (size_t i=0; i<length; ++i)
    {
        size_t index = ( n_indices == 0 ? i : indices[i]);
        v[index] -= storedSlidingFactor;
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void VectorSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode< RbVector<double> >* >(newN) ;
    
}


void VectorSlideProposal::setProposalTuningParameter(double tp)
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
void VectorSlideProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        lambda /= (2.0 - rate/0.44 );
    }
    
}

