#include "ElementScaleProposal.h"
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
ElementScaleProposal::ElementScaleProposal( std::vector< StochasticNode<double> *> n, double l) : Proposal(),
    variables( n ),
    lambda( l ),
    length( variables.size() ),
    storedScalingFactor( 1.0 )
{
    // tell the base class to add the node
    for (std::vector< StochasticNode<double> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
    {
        addNode( *it );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ElementScaleProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ElementScaleProposal* ElementScaleProposal::clone( void ) const
{
    
    return new ElementScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& ElementScaleProposal::getProposalName( void ) const
{
    static std::string name = "VectorScaling";
    
    return name;
}


double ElementScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double ElementScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    
    // copy value
    storedScalingFactor = scaling_factor;
    
    for (size_t index=0; index<length; ++index)
    {
        variables[index]->getValue() *= scaling_factor;
    }
    
    // compute the Hastings ratio
    double lnHastingsratio = length * log( scaling_factor );
    
    return lnHastingsratio;
}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void ElementScaleProposal::prepareProposal( void )
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
void ElementScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void ElementScaleProposal::undoProposal( void )
{
    
    for (size_t index=0; index<length; ++index)
    {
        variables[index]->getValue() /= storedScalingFactor;
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void ElementScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
}


void ElementScaleProposal::setProposalTuningParameter(double tp)
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
void ElementScaleProposal::tune( double rate )
{
    
    if ( rate > 0.23 )
    {
        lambda *= (1.0 + ((rate-0.23)/0.77) );
    }
    else
    {
        lambda /= (2.0 - rate/0.23 );
    }
    
}

