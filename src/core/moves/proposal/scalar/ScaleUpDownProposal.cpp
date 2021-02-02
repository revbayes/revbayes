#include <cmath>
#include <iostream>

#include "DistributionNormal.h"
#include "ScaleUpDownProposal.h"
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
ScaleUpDownProposal::ScaleUpDownProposal( ContinuousStochasticNode *n1, ContinuousStochasticNode *n2, double l) : Proposal(),
    variable1( n1 ),
    variable2( n2 ),
    storedValue1( 0.0 ),
    storedValue2( 0.0 ),
    lambda( l )
{
    // tell the base class to add the node
    addNode( variable1 );
    addNode( variable2 );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ScaleUpDownProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ScaleUpDownProposal* ScaleUpDownProposal::clone( void ) const
{
    
    return new ScaleUpDownProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& ScaleUpDownProposal::getProposalName( void ) const
{
    static std::string name = "ScaleUpDown";
    
    return name;
}


double ScaleUpDownProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and ScaleUpDowns the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double ScaleUpDownProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    double sv1 = variable1->getValue();
    double sv2 = variable2->getValue();
    
    storedValue1 = sv1;
    storedValue2 = sv2;
    
    double u = rng->uniform01();
    double c = exp( lambda * (u - 0.5) );
    
    c = 1.0;
    
    sv1 *= c;
    sv2 /= c;
    
    variable1->setValue( new double(sv1) );
    variable2->setValue( new double(sv2) );
    
    //	double pr = (1.0 - 1.0 - 2.0) * log(c);
    return 0.0;
}


/**
 *
 */
void ScaleUpDownProposal::prepareProposal( void )
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
void ScaleUpDownProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void ScaleUpDownProposal::undoProposal( void )
{
    // swap current value and stored value
    variable1->setValue( new double(storedValue1) );
    variable2->setValue( new double(storedValue2) );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void ScaleUpDownProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == variable1)
    {
        variable1 = static_cast<ContinuousStochasticNode* >( newN );
    }
    
    if (oldN == variable2)
    {
        variable2 = static_cast<ContinuousStochasticNode* >( newN );
    }
    
}


void ScaleUpDownProposal::setProposalTuningParameter(double tp)
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
void ScaleUpDownProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        lambda /= (2.0 - rate/0.44 );
    }
    
    lambda = fmin(10000, lambda);
    
}

