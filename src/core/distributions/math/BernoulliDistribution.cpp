#include "BernoulliDistribution.h"

#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/*
 * Bernoulli Distribution Constructor
 *
 * @param p a value of type double for the probability of succcess
 *
 */

BernoulliDistribution::BernoulliDistribution(const TypedDagNode<double> *p) : TypedDistribution<std::int64_t>( new std::int64_t( 0 ) ),
    p( p )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( p );
    
    double u = GLOBAL_RNG->uniform01();
    *value = u > p->getValue() ? 0 : 1;
}


BernoulliDistribution::~BernoulliDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BernoulliDistribution* BernoulliDistribution::clone( void ) const {
    
    return new BernoulliDistribution( *this );
}


double BernoulliDistribution::computeLnProbability( void )
{

    // check that the value is actually a bernoulli trial (1 or 0)
    if ( *value != 1 && *value != 0 )
    {
        return RbConstants::Double::neginf;
    }
    
    if ( *value > 1 || *value < 0 )
        return RbConstants::Double::neginf;

    return log( *value == 0 ? ( 1 - p->getValue() ) : p->getValue() );
}


void BernoulliDistribution::redrawValue( void ) {
    
    double u = GLOBAL_RNG->uniform01();
    *value = u > p->getValue() ? 0 : 1;
    
}


/** Swap a parameter of the distribution */
void BernoulliDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == p) 
    {
        p = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}


