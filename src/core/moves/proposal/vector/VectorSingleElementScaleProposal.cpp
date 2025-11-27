#include "VectorSingleElementScaleProposal.h"

#include <cmath>
#include <ostream>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


VectorSingleElementScaleProposal::VectorSingleElementScaleProposal( StochasticNode< RbVector<double> > *n, double l) : Proposal(),
    variable( n ),
    lambda( l ),
    index( 0 ),
    storedValue( 0.0 )
{
    // tell the base class to add the node
    addNode( variable );
}


void VectorSingleElementScaleProposal::cleanProposal( void )
{
    variable->clearTouchedElementIndices();
}


VectorSingleElementScaleProposal* VectorSingleElementScaleProposal::clone( void ) const
{
    return new VectorSingleElementScaleProposal( *this );
}


const std::string& VectorSingleElementScaleProposal::getProposalName( void ) const
{
    static std::string name = "ElementScaling";
    return name;
}


double VectorSingleElementScaleProposal::getProposalTuningParameter( void ) const
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
double VectorSingleElementScaleProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    RbVector<double> &val = variable->getValue();

    // choose an index
    index = size_t(rng->uniform01() * val.size());

    // copy value
    storedValue = val[index];

    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    val[index] *= scaling_factor;

    variable->addTouchedElementIndex(index);

    // compute the Hastings ratio
    double lnHastingsratio = log( scaling_factor );

    return lnHastingsratio;
}


void VectorSingleElementScaleProposal::prepareProposal( void ) {}


void VectorSingleElementScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }
}


void VectorSingleElementScaleProposal::undoProposal( void )
{
    std::vector<double>& v = variable->getValue();
    v[index] = storedValue;
    variable->clearTouchedElementIndices();
}


void VectorSingleElementScaleProposal::swapNodeInternal(DagNode *, DagNode *newN)
{
    variable = static_cast<StochasticNode< RbVector<double> >* >(newN) ;
}


void VectorSingleElementScaleProposal::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.23.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 *
 * @param rate current acceptance rate
 */
void VectorSingleElementScaleProposal::tune( double rate )
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

