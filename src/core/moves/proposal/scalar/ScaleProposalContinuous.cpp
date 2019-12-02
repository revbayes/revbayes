#include "ScaleProposalContinuous.h"

#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "ContinuousStochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** 
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
ScaleProposalContinuous::ScaleProposalContinuous( ContinuousStochasticNode *n, double l, double p) : Proposal(p),
    variable( n ),
    storedValue( 0.0 ),
    lambda( l )
{
    // tell the base class to add the node
    addNode( variable );
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ScaleProposalContinuous::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal. 
 */
ScaleProposalContinuous* ScaleProposalContinuous::clone( void ) const 
{
    
    return new ScaleProposalContinuous( *this );
}


/**
 * Get Proposals' name of object 
 *
 * \return The Proposals' name.
 */
const std::string& ScaleProposalContinuous::getProposalName( void ) const 
{
    static std::string name = "Scaling";
    
    return name;
}


double ScaleProposalContinuous::getProposalTuningParameter( void ) const
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
double ScaleProposalContinuous::doProposal( void )
{
    
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    double &val = variable->getValue();
        
    // copy value
    storedValue = val;
    
//    double min = variable->getMin();
//    double max = variable->getMax();
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    double newVal = val * scaling_factor;
    
    /* reflect the new value */
//    while ( newVal < min || newVal > max ) {
//        if ( newVal < min )
//        {
//            newVal = pow(min, 2) / newVal;
//        }
//        else if ( newVal > max )
//        {
//            newVal = pow(max, 2) / newVal;
//        }
//    }
    
    // compute the Hastings ratio
    double lnHastingsratio = log( newVal / val ); // Jacobian
    val = newVal;
    
    return lnHastingsratio;
}


/**
 *
 */
void ScaleProposalContinuous::prepareProposal( void )
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
void ScaleProposalContinuous::printParameterSummary(std::ostream &o, bool name_only) const
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
void ScaleProposalContinuous::undoProposal( void ) 
{
    // swap current value and stored value
    variable->setValue( new double(storedValue) );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void ScaleProposalContinuous::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< ContinuousStochasticNode* >(newN) ;
    
}


void ScaleProposalContinuous::setProposalTuningParameter(double tp)
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
void ScaleProposalContinuous::tune( double rate ) 
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else 
    {
        lambda /= (2.0 - rate/p);
    }
    
}

