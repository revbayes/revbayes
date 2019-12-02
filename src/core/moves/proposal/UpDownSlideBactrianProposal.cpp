#include "DistributionNormal.h"
#include "UpDownSlideBactrianProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
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
UpDownSlideBactrianProposal::UpDownSlideBactrianProposal( double l ) : Proposal(),
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
void UpDownSlideBactrianProposal::addVariable( StochasticNode<double> *v, bool up)
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
void UpDownSlideBactrianProposal::addVariable( StochasticNode<RbVector<double> > *v, bool up)
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
void UpDownSlideBactrianProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
UpDownSlideBactrianProposal* UpDownSlideBactrianProposal::clone( void ) const
{
    
    return new UpDownSlideBactrianProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& UpDownSlideBactrianProposal::getProposalName( void ) const
{
    static std::string name = "UpDownSlideBactrian";
    
    return name;
}


double UpDownSlideBactrianProposal::getProposalTuningParameter( void ) const
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
double UpDownSlideBactrianProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
        
    double u = rng->uniform01();
    double delta = RbStatistics::Normal::rv(*GLOBAL_RNG) * RbConstants::BACT_SD;
    
    // See Yang and Rodriguez (2013) SI eqns 19 and 20
    // Currently hard-coding m = 0.95
    if ( u < 0.5 ) {
        delta += 0.95;
    } else {
        delta -= 0.95;
    }
    
    delta = ( lambda * delta );
    
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
    
    // this is a symmetric proposal so the hasting ratio is 0.0
    return 0.0;

}


/**
 * Prepare the proposal, e.g., pick the UpDown that we want to change.
 * Here we do not need to do any preparation.
 */
void UpDownSlideBactrianProposal::prepareProposal( void )
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
void UpDownSlideBactrianProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void UpDownSlideBactrianProposal::removeVariable( StochasticNode<double> *v, bool up)
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
void UpDownSlideBactrianProposal::removeVariable( StochasticNode<RbVector<double> > *v, bool up)
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
void UpDownSlideBactrianProposal::undoProposal( void )
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
void UpDownSlideBactrianProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
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


void UpDownSlideBactrianProposal::setProposalTuningParameter(double tp)
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
void UpDownSlideBactrianProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        lambda /= (2.0 - rate/0.44 );
    }
    
    if ( lambda > 1 ) {
        lambda = fmin(10000, lambda);
    } else {
        lambda = fmax(1/10000, lambda);
    }
    
}

