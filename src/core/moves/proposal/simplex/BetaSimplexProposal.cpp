#include "BetaSimplexProposal.h"

#include <cstddef>
#include <cmath>
#include <iostream>

#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbConstants.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
BetaSimplexProposal::BetaSimplexProposal( StochasticNode<Simplex> *n, double a, double p ) : SimpleProposal<Simplex>( n, p ),
    storedValue( Simplex() ),
    alpha( a )
{
    // tell the base class to add the node
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void BetaSimplexProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
BetaSimplexProposal* BetaSimplexProposal::clone( void ) const
{
    
    return new BetaSimplexProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& BetaSimplexProposal::getProposalName( void ) const
{
    static std::string name = "BetaSimplex";
    
    return name;
}


double BetaSimplexProposal::getProposalTuningParameter( void ) const
{
    return alpha;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new simplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double BetaSimplexProposal::propose( Simplex &value )
{
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // we need to know the number of categories
    size_t cats = value.size();
    
    // store the value
    storedValue = value;
    
    // randomly draw a new index
    size_t chosen_index = size_t( floor(rng->uniform01()*double(cats)) );
    double current_value = value[chosen_index];
    
    // draw new rates and compute the hastings ratio at the same time
    double a = alpha + 1.0;
    double b = (a - 1.0) / current_value - a + 2.0;
    double new_value = RbStatistics::Beta::rv(a, b, *rng);
    
    // set the value
    value[chosen_index] = new_value;

    double scaling_factor_other_values = (1.0 - new_value) / (1.0 - current_value);
    
    double sum = 0.0;
    double ln_Hastings_ratio = 0.0;
    for (size_t i = 0; i < cats; i++)
    {
        if ( i != chosen_index )
        {
            value[i] *= scaling_factor_other_values;
        }
        
        // test for 0-values
        if ( value[i] < 1E-100)
        {
            return RbConstants::Double::neginf;
        }
        
        sum += value[i];
    }
    
    // mrm 08-07-2017: I've noticed the sum of this
    // simplex can gradually drift until it the sum
    // is substantially different from one.
    // normalize the vector to enforce sum = 1
    for (size_t i = 0; i < cats; ++i)
    {
        value[i] /= sum;
    }
    
    try
    {
        // compute the Hastings ratio
        double forward = RbStatistics::Beta::lnPdf(a, b, new_value);
        double new_a = alpha + 1.0;
        double new_b = (a-1.0) / value[chosen_index] - a + 2.0;
        double backward = RbStatistics::Beta::lnPdf(new_a, new_b, current_value);
        
        ln_Hastings_ratio = backward - forward;
        
        // include the Jacobian for the scaling of the other values
        ln_Hastings_ratio += (cats - 2) * log(scaling_factor_other_values) - (cats - 1) * log(sum);
    }
	catch (RbException &e)
    {
        ln_Hastings_ratio = RbConstants::Double::neginf;
    }

    return ln_Hastings_ratio;
}


/**
 *
 */
void BetaSimplexProposal::prepareProposal( void )
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
void BetaSimplexProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void BetaSimplexProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->setValue( storedValue.clone() );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void BetaSimplexProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Simplex>* >(newN) ;
    
}


void BetaSimplexProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void BetaSimplexProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    
    if ( rate > p )
    {
        alpha /= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        alpha *= (2.0 - rate/p);
    }
    
    // set a hard maximum of 100
    alpha = fmin(100, alpha);
}

