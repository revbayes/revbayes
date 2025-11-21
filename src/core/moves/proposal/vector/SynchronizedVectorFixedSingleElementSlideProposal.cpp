#include "SynchronizedVectorFixedSingleElementSlideProposal.h"

#include <cstddef>

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
SynchronizedVectorFixedSingleElementSlideProposal::SynchronizedVectorFixedSingleElementSlideProposal( const std::vector< StochasticNode< RbVector<double> > *> &n, double l, const std::vector<std::int64_t> &i) : Proposal(),
    variables( n ),
    lambda( l ),
    indices( i ),
    stored_delta( 0.0 )
{
    // tell the base class to add the node
    for (size_t i=0; i<variables.size(); ++i)
    {
        addNode( variables[i] );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SynchronizedVectorFixedSingleElementSlideProposal::cleanProposal( void )
{

    for (size_t i=0; i<variables.size(); ++i)
    {
        variables[i]->clearTouchedElementIndices();
    }

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SynchronizedVectorFixedSingleElementSlideProposal* SynchronizedVectorFixedSingleElementSlideProposal::clone( void ) const
{
    
    return new SynchronizedVectorFixedSingleElementSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SynchronizedVectorFixedSingleElementSlideProposal::getProposalName( void ) const
{
    static std::string name = "SynchronizedVectorFixedSingleElementSliding";
    
    return name;
}


double SynchronizedVectorFixedSingleElementSlideProposal::getProposalTuningParameter( void ) const
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
double SynchronizedVectorFixedSingleElementSlideProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double delta  = ( lambda * ( u - 0.5 ) );
    stored_delta = delta;
    
    for (size_t i=0; i<variables.size(); ++i)
    {
        RbVector<double> &val = variables[i]->getValue();
        for (size_t j=0; j<indices.size(); ++j)
        {
            size_t index = indices[j];
            val[index] += delta;
            variables[i]->addTouchedElementIndex(index);
        }
    }
    
    
    return 0.0;
}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void SynchronizedVectorFixedSingleElementSlideProposal::prepareProposal( void )
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
void SynchronizedVectorFixedSingleElementSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void SynchronizedVectorFixedSingleElementSlideProposal::undoProposal( void )
{
    
    for (size_t i=0; i<variables.size(); ++i)
    {
        RbVector<double> &val = variables[i]->getValue();
        for (size_t j=0; j<indices.size(); ++j)
        {
            size_t index = indices[j];
            val[index] -= stored_delta;
            variables[i]->clearTouchedElementIndices();
        }
        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SynchronizedVectorFixedSingleElementSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<RbVector<double> > *>(newN);
        }
    }
    
}


void SynchronizedVectorFixedSingleElementSlideProposal::setProposalTuningParameter(double tp)
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
void SynchronizedVectorFixedSingleElementSlideProposal::tune( double rate )
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

