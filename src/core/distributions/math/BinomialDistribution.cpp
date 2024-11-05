#include "BinomialDistribution.h"

#include "DistributionBinomial.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** BinomialDistribution Constructor
 * @param n number of trials
 * @param p the probability
*/
BinomialDistribution::BinomialDistribution(const TypedDagNode<std::int64_t> *n, const TypedDagNode<double> *p) : TypedDistribution<std::int64_t>( new long( 0 ) ),
    n( n ),
    p( p )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( n );
    addParameter( p );
    
    *value = RbStatistics::Binomial::rv(n->getValue(), p->getValue(), *GLOBAL_RNG);
}


BinomialDistribution::~BinomialDistribution(void)
{

    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BinomialDistribution* BinomialDistribution::clone( void ) const
{

    return new BinomialDistribution( *this );
}


double BinomialDistribution::computeLnProbability( void )
{
    
    // check that the value is inside the boundaries
    if ( *value > n->getValue() || *value < 0 || n->getValue() < 0)
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Binomial::lnPdf(n->getValue(), p->getValue(), *value);
}



void BinomialDistribution::redrawValue( void )
{

    *value = RbStatistics::Binomial::rv(n->getValue(), p->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void BinomialDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == p)
    {
        p = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == n)
    {
        n = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }

}


