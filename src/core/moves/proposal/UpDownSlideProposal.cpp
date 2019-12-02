#include "UpDownSlideProposal.h"
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
UpDownSlideProposal::UpDownSlideProposal( double l ) : Proposal(),
    upScalarVariables(  ),
    upVectorVariables(  ),
    downScalarVariables(  ),
    downVectorVariables(  ),
    lambda( l ),
    storedDelta( 0.0 )
{
    
}


/**
 * Add a variable.
 */
void UpDownSlideProposal::addVariable( StochasticNode<double> *v, bool up)
{
    
    
    // add it to the nodes vector
    addNode( v );
    
    // add it to my internal vector
    if ( up == true )
    {
        upScalarVariables.push_back( v );
    }
    else
    {
        downScalarVariables.push_back( v );
    }
    
}


/**
 * Add a variable.
 */
void UpDownSlideProposal::addVariable( StochasticNode<RbVector<double> > *v, bool up)
{
    
    
    // add it to the nodes vector
    addNode( v );
    
    // add it to my internal vector
    if ( up == true )
    {
        upVectorVariables.push_back( v );
    }
    else
    {
        downVectorVariables.push_back( v );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void UpDownSlideProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
UpDownSlideProposal* UpDownSlideProposal::clone( void ) const
{
    
    return new UpDownSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& UpDownSlideProposal::getProposalName( void ) const
{
    static std::string name = "UpDownSlide";
    
    return name;
}


double UpDownSlideProposal::getProposalTuningParameter( void ) const
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
double UpDownSlideProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double delta = lambda * ( u - 0.5 ) ;
    
    // copy value
    storedDelta = delta;
    
    // Slide all the single variable up
    for (size_t i=0; i<upScalarVariables.size(); ++i)
    {
        upScalarVariables[i]->getValue() += delta;
    }
    
    // Slide all the vector variables up
    for (size_t i=0; i<upVectorVariables.size(); ++i)
    {
        RbVector<double> &v = upVectorVariables[i]->getValue();
        for (size_t j=0; j<v.size(); ++j)
        {
            v[j] += delta;
        }
    }
    
    // Slide all the single variable down
    for (size_t i=0; i<downScalarVariables.size(); ++i)
    {
        downScalarVariables[i]->getValue() -= delta;
    }
    
    // Slide all the vector variables down
    for (size_t i=0; i<downVectorVariables.size(); ++i)
    {
        RbVector<double> &v = downVectorVariables[i]->getValue();
        for (size_t j=0; j<v.size(); ++j)
        {
            v[j] -= delta;
        }

    }
    
    // compute the Hastings ratio (symmetric move)
    double lnHastingsratio = 0.0;
    
    return lnHastingsratio;
}


/**
 * Prepare the proposal, e.g., pick the UpDown that we want to change.
 * Here we do not need to do any preparation.
 */
void UpDownSlideProposal::prepareProposal( void )
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
void UpDownSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }
    
}


/**
 * Remove a variable.
 */
void UpDownSlideProposal::removeVariable( StochasticNode<double> *v, bool up)
{
    
    
    // add it to the nodes vector
    removeNode( v );
    
    // add it to my internal vector
    if ( up == true )
    {
        upScalarVariables.erase(std::remove(upScalarVariables.begin(), upScalarVariables.end(), v), upScalarVariables.end());
    }
    else
    {
        downScalarVariables.erase(std::remove(downScalarVariables.begin(), downScalarVariables.end(), v), downScalarVariables.end());
    }
    
}


/**
 * Remove a variable.
 */
void UpDownSlideProposal::removeVariable( StochasticNode<RbVector<double> > *v, bool up)
{
    
    
    // add it to the nodes vector
    removeNode( v );
    
    // add it to my internal vector
    if ( up == true )
    {
        upVectorVariables.erase(std::remove(upVectorVariables.begin(), upVectorVariables.end(), v), upVectorVariables.end());
    }
    else
    {
        downVectorVariables.erase(std::remove(downVectorVariables.begin(), downVectorVariables.end(), v), downVectorVariables.end());
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void UpDownSlideProposal::undoProposal( void )
{
    
    double delta = -storedDelta;
    
    // Slide all the single variable up
    for (size_t i=0; i<upScalarVariables.size(); ++i)
    {
        upScalarVariables[i]->getValue() += delta;
    }
    
    // Slide all the vector variables up
    for (size_t i=0; i<upVectorVariables.size(); ++i)
    {
        RbVector<double> &v = upVectorVariables[i]->getValue();
        for (size_t j=0; j<v.size(); ++j)
        {
            v[j] += delta;
        }
    }
    
    // Slide all the single variable down
    for (size_t i=0; i<downScalarVariables.size(); ++i)
    {
        downScalarVariables[i]->getValue() -= delta;
    }
    
    // Slide all the vector variables down
    for (size_t i=0; i<downVectorVariables.size(); ++i)
    {
        RbVector<double> &v = downVectorVariables[i]->getValue();
        for (size_t j=0; j<v.size(); ++j)
        {
            v[j] -= delta;
        }
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void UpDownSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < upScalarVariables.size(); ++i)
    {
        if ( upScalarVariables[i] == oldN )
        {
            upScalarVariables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
    for (size_t i = 0; i < upVectorVariables.size(); ++i)
    {
        if ( upVectorVariables[i] == oldN )
        {
            upVectorVariables[i] = static_cast<StochasticNode<RbVector<double> > *>(newN);
        }
    }
    
    for (size_t i = 0; i < downScalarVariables.size(); ++i)
    {
        if ( downScalarVariables[i] == oldN )
        {
            downScalarVariables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
    for (size_t i = 0; i < downVectorVariables.size(); ++i)
    {
        if ( downVectorVariables[i] == oldN )
        {
            downVectorVariables[i] = static_cast<StochasticNode<RbVector<double> > *>(newN);
        }
    }
    
}


void UpDownSlideProposal::setProposalTuningParameter(double tp)
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
void UpDownSlideProposal::tune( double rate )
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

