#include "MultinomialDistribution.h"

#include <cstddef>
#include <cstdint>

#include "DistributionMultinomial.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*Multinomial Distribution Constructor
 * @param p A simplex of the the probabilities for each category
 * @param n A std::int64_t for the number of trials
 */

MultinomialDistribution::MultinomialDistribution(const TypedDagNode<Simplex> *p, const TypedDagNode<std::int64_t> *n) : TypedDistribution< RbVector<std::int64_t> >( new RbVector<std::int64_t>() ),
    p( p ),
    n( n )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( p );
    addParameter( n );
        
    *value = RbStatistics::Multinomial::rv(p->getValue(), size_t( n->getValue() ), *GLOBAL_RNG);
}


MultinomialDistribution::~MultinomialDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
    // The model needs to do that!
}



MultinomialDistribution* MultinomialDistribution::clone( void ) const
{
    return new MultinomialDistribution( *this );
}


double MultinomialDistribution::computeLnProbability( void )
{
    return RbStatistics::Multinomial::lnPdf(p->getValue(), *value);
}


void MultinomialDistribution::redrawValue( void )
{
    *value = RbStatistics::Multinomial::rv(p->getValue(), size_t( n->getValue() ), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void MultinomialDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == p)
    {
        p = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == n)
    {
        n = static_cast<const TypedDagNode<std::int64_t>* >(newP);
    }
}
