#include "VectorScaleProposal.h"

#include <cmath>
#include <ostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


VectorScaleProposal::VectorScaleProposal( StochasticNode< RbVector<double> > *n, double l) : Proposal(),
    variable( n ),
    lambda( l ),
    length( variable->getValue().size() ),
    storedScalingFactor( 1.0 )
{
    // tell the base class to add the node
    addNode( variable );    
}


void VectorScaleProposal::cleanProposal( void )
{
    variable->clearTouchedElementIndices();
}


VectorScaleProposal* VectorScaleProposal::clone( void ) const
{   
    return new VectorScaleProposal( *this );
}


const std::string& VectorScaleProposal::getProposalName( void ) const
{
    static std::string name = "VectorScaling";    
    return name;
}


double VectorScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A scaling proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the proposal.
 *
 * \return The log hastings ratio.
 */
double VectorScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    RbVector<double> &val = variable->getValue();
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    
    // copy value
    storedScalingFactor = scaling_factor;
    
    for (size_t index=0; index<length; ++index)
    {
        val[index] *= scaling_factor;
    }
    
    // compute the Hastings ratio
    double lnHastingsratio = length * log( scaling_factor );
    
    return lnHastingsratio;
}


void VectorScaleProposal::prepareProposal( void ) {}


void VectorScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{    
    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }   
}



void VectorScaleProposal::undoProposal( void )
{
    RbVector<double>& v = variable->getValue();
    
    for (size_t index=0; index<length; ++index)
    {
        v[index] /= storedScalingFactor;
    }    
}


void VectorScaleProposal::swapNodeInternal(DagNode *, DagNode *newN)
{   
    variable = static_cast<StochasticNode< RbVector<double> >* >(newN) ;   
}


void VectorScaleProposal::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 *
 * @param rate current acceptance rate
 */
void VectorScaleProposal::tune( double rate )
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

