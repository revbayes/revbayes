#include <cmath>
#include <iostream>

#include "DistributionNormal.h"
#include "LevyJumpProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "ContinuousStochasticNode.h"
#include "Proposal.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
LevyJumpProposal::LevyJumpProposal( ContinuousStochasticNode *n, double d) : Proposal(),
    variable( n ),
    storedValue( 0.0 ),
    delta( d )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void LevyJumpProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
LevyJumpProposal* LevyJumpProposal::clone( void ) const
{
    
    return new LevyJumpProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& LevyJumpProposal::getProposalName( void ) const
{
    static std::string name = "Sliding";
    
    return name;
}


double LevyJumpProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double LevyJumpProposal::doProposal( void )
{
    
    // Gottardo and Rafferty 2009
    
    
    // store variable for restoration
    double sv = variable->getValue();
    storedValue = sv;
    
    // turn on/off variable
    double u = GLOBAL_RNG->uniform01();
    double p = 0.5;
    
    // propose new value
    double nv = RbStatistics::Normal::rv(sv, delta, *GLOBAL_RNG);
    if (u < p)
        nv = 0.0;
    
    variable->setValue( new double(nv) );
    
    // MH ratio
    double probFwd = p * (nv==0.0) + (1-p) * RbStatistics::Normal::pdf(sv, delta, nv) * (nv!=0.0);
    double probRev = p * (sv==0.0) + (1-p) * RbStatistics::Normal::pdf(nv, delta, sv) * (sv!=0.0);
    
    return log( probRev / probFwd );
}


/**
 *
 */
void LevyJumpProposal::prepareProposal( void )
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
void LevyJumpProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void LevyJumpProposal::undoProposal( void )
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
void LevyJumpProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< ContinuousStochasticNode* >(newN) ;
    
}


void LevyJumpProposal::setProposalTuningParameter(double tp)
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
void LevyJumpProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        delta *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        delta /= (2.0 - rate/0.44 );
    }
    
}

