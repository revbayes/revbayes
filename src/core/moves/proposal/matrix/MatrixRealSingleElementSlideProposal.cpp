#include "MatrixRealSingleElementSlideProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MatrixRealSingleElementSlideProposal::MatrixRealSingleElementSlideProposal( StochasticNode<MatrixReal> *n, double l, bool s) : Proposal(),
    array(NULL),
    matrix( n ),
    delta( l ),
    indexa(0),
    indexb(0),
    storedValue( 0.0 ),
    symmetric( s )
{
    // tell the base class to add the node
    addNode( matrix );

}


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MatrixRealSingleElementSlideProposal::MatrixRealSingleElementSlideProposal( StochasticNode<RbVector<RbVector<double> > > *n, double l, bool s) : Proposal(),
    array( n ),
    matrix(NULL),
    delta( l ),
    indexa(0),
    indexb(0),
    storedValue( 0.0 ),
    symmetric( s )
{
    // tell the base class to add the node
    addNode( array );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MatrixRealSingleElementSlideProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MatrixRealSingleElementSlideProposal* MatrixRealSingleElementSlideProposal::clone( void ) const
{
    
    return new MatrixRealSingleElementSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MatrixRealSingleElementSlideProposal::getProposalName( void ) const
{
    static std::string name = "MatrixRealSingleElementSlideMove";
    
    return name;
}


double MatrixRealSingleElementSlideProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double MatrixRealSingleElementSlideProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scalingFactor = delta * ( u - 0.5 );
    
    if (array != NULL)
    {
        RbVector<RbVector<double> >& v = array->getValue();
        // choose an index
        indexa = size_t( rng->uniform01() * v.size() );
        indexb = size_t( rng->uniform01() * v.front().size() );

        // copy value
        storedValue = v[indexa][indexb];

        v[indexa][indexb] += scalingFactor;

        if ( symmetric == true && indexa != indexb)
        {
            v[indexb][indexa] += scalingFactor;
        }

        array->addTouchedElementIndex(indexa*v.size() + indexb);
    }
    else
    {
        MatrixReal& v = matrix->getValue();
        // choose an index
        indexa = size_t( rng->uniform01() * v.getNumberOfRows() );
        indexb = size_t( rng->uniform01() * v.getNumberOfColumns() );

        // copy value
        storedValue = v[indexa][indexb];

        v[indexa][indexb] += scalingFactor;

        if ( symmetric == true && indexa != indexb)
        {
            v[indexb][indexa] += scalingFactor;
        }

        matrix->addTouchedElementIndex(indexa*v.getNumberOfRows() + indexb);
    }
    
    // this is a symmetric proposal so the hasting ratio is 0.0
    return 0.0;
}


/**
 *
 */
void MatrixRealSingleElementSlideProposal::prepareProposal( void )
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
void MatrixRealSingleElementSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = ";
    if (name_only == false)
    {
        o << delta;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MatrixRealSingleElementSlideProposal::undoProposal( void )
{
    if (array != NULL)
    {
        RbVector<RbVector<double> >& v = array->getValue();
        v[indexa][indexb] = storedValue;
        if ( symmetric == true )
        {
            v[indexb][indexa] = storedValue;
        }
        array->clearTouchedElementIndices();
    }
    else
    {
        MatrixReal& v = matrix->getValue();
        v[indexa][indexb] = storedValue;
        if ( symmetric == true )
        {
            v[indexb][indexa] = storedValue;
        }
        matrix->clearTouchedElementIndices();
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MatrixRealSingleElementSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if (oldN == array)
    {
        array = static_cast< StochasticNode<RbVector<RbVector<double> > >* >(newN) ;
    }
    else
    {
        matrix = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    }
    
}


void MatrixRealSingleElementSlideProposal::setProposalTuningParameter(double tp)
{
    delta = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void MatrixRealSingleElementSlideProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        delta *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        delta /= (2.0 - rate/0.44 );
    }
    
    delta = fmin(10000, delta);
    
}

