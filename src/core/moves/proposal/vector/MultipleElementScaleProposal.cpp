#include <cmath>
#include <iostream>
#include <cstddef>
#include <vector>

#include "DistributionNormal.h"
#include "MultipleElementScaleProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbStatisticsHelper.h"
#include "RbException.h"
#include "Proposal.h"
#include "RbConstants.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MultipleElementScaleProposal::MultipleElementScaleProposal( std::vector< StochasticNode<double> *> n, size_t k, double l) : Proposal(),
    variables( n ),
    numToMove( k ),
    lambda( l ),
    length( variables.size() ),
    storedScalingFactors( std::vector<double>(numToMove,1.0) ),
    storedChosenElements( std::vector<size_t>(numToMove,0) )
{
    // tell the base class to add the node
    for (std::vector< StochasticNode<double> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
    {
        addNode( *it );
    }
    
    if (numToMove < 1)
    {
        throw(RbException("MultipleElementScale must move at least one vector element at a time"));
    }

    if (numToMove > length)
    {
        throw(RbException("MultipleElementScale cannot move more vector elements than there are elements in the vector"));
    }
        
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MultipleElementScaleProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MultipleElementScaleProposal* MultipleElementScaleProposal::clone( void ) const
{
    
    return new MultipleElementScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MultipleElementScaleProposal::getProposalName( void ) const
{
    static std::string name = "MultipleVectorScaling";
    
    return name;
}


double MultipleElementScaleProposal::getProposalTuningParameter( void ) const
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
double MultipleElementScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // Choose which elements to move
    std::vector<size_t> to_move;
    std::vector<size_t> tmpV;
    
    for (size_t i=0; i<length; i++)
    {
        tmpV.push_back(i);
    }

    RbStatistics::Helper::randomlySelectFromVectorWithoutReplacement<size_t>(tmpV, to_move, numToMove, *rng);
    
    double lnHastingsratio = 0.0;
    
    size_t index;
    
    // Generate new values (no reflection, so we simply abort later if we propose value here outside of support)
    for (size_t i=0; i<numToMove; ++i)
    {
        // Bactrian RV
        double u = rng->uniform01();
        double delta = RbStatistics::Normal::rv(*GLOBAL_RNG) * RbConstants::BACT_SD;
        
        // See Yang and Rodriguez (2013) SI eqns 19 and 20
        // Currently hard-coding m = 0.95
        if ( u < 0.5 ) {
            delta += 0.95;
        } else {
            delta -= 0.95;
        }
        
        double scaling_factor = lambda * delta;
        lnHastingsratio += scaling_factor;
        scaling_factor = std::exp( scaling_factor );
        
        index = to_move[i];
        
        // copy value
        storedScalingFactors[i] = scaling_factor;
        storedChosenElements[i] = index;
        
        variables[index]->getValue() *= scaling_factor;
        
    }
        
    return lnHastingsratio;
}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void MultipleElementScaleProposal::prepareProposal( void )
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
void MultipleElementScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void MultipleElementScaleProposal::undoProposal( void )
{
    
    for (size_t i=0; i<numToMove; ++i)
    {
        size_t index = storedChosenElements[i];
        variables[index]->getValue() /= storedScalingFactors[i];
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MultipleElementScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
}


void MultipleElementScaleProposal::setProposalTuningParameter(double tp)
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
void MultipleElementScaleProposal::tune( double rate )
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

