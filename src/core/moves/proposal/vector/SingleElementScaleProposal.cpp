#include "SingleElementScaleProposal.h"

#include <cmath>
#include <ostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

SingleElementScaleProposal::SingleElementScaleProposal( std::vector<StochasticNode<double> *> n, double l) : Proposal(),
    variables( n ),
    lambda( l ),
    index( 0 ),
    storedValue( 0.0 )
{
    // tell the base class to add the node
    for (std::vector< StochasticNode<double> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
    {
        addNode( *it );
    }
    
}


SingleElementScaleProposal* SingleElementScaleProposal::clone( void ) const
{    
    return new SingleElementScaleProposal( *this );
}


const std::string& SingleElementScaleProposal::getProposalName( void ) const
{
    static std::string name = "ElementScaling";
    
    return name;
}


double SingleElementScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A scaling proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current value by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the proposal.
 *
 * \return The log hastings ratio.
 */
double SingleElementScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // choose an index
    index = size_t(rng->uniform01() * variables.size());
    
    double &val = variables[index]->getValue();
    
    // copy value
    storedValue = val;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    val *= scaling_factor;
    
    // compute the Hastings ratio
    double lnHastingsratio = log( scaling_factor );
    
    return lnHastingsratio;
}


void SingleElementScaleProposal::prepareProposal( void ) {}

void SingleElementScaleProposal::cleanProposal( void ) {}


void SingleElementScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }
    
}


void SingleElementScaleProposal::undoProposal( void )
{
    double &v = variables[index]->getValue();
    v = storedValue;    
}


void SingleElementScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
}


void SingleElementScaleProposal::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 *
 * @param rate current acceptance rate
 */
void SingleElementScaleProposal::tune( double rate )
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

