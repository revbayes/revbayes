#include <vector>
#include <cstdint>

#include "BivariatePoissonDistribution.h"
#include "DistributionBivariatePoisson.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "Cloner.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

BivariatePoissonDistribution::BivariatePoissonDistribution(const TypedDagNode< double > *theta1x, const TypedDagNode< double > *theta2x, const TypedDagNode< double > *theta0x) : TypedDistribution< RbVector<std::int64_t> >( new RbVector<std::int64_t>() ),
  theta1( theta1x ),
  theta2( theta2x ),
  theta0( theta0x )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
  addParameter( theta1 );
  addParameter( theta2 );
  addParameter( theta0 );

  //*value = RbStatistics::BivariatePoisson::rv(theta1->getValue(), theta2->getValue(), theta0->getValue(), *GLOBAL_RNG);
  redrawValue();
}


BivariatePoissonDistribution::~BivariatePoissonDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


BivariatePoissonDistribution* BivariatePoissonDistribution::clone( void ) const {
    return new BivariatePoissonDistribution( *this );
}


double BivariatePoissonDistribution::computeLnProbability( void )
{

  std::vector<std::int64_t> v = *value;
  std::int64_t x1, x2;
  x1 = v.front();
  x2 = v.front();

    // check that the value is inside the boundaries
    if ( x1 < 0 || x2 < 0 )
    {
        return RbConstants::Double::neginf;
    }

    return RbStatistics::BivariatePoisson::lnPdf(theta1->getValue(), theta2->getValue(), theta0->getValue(), *value);
}


void BivariatePoissonDistribution::redrawValue( void ) {
  *value = RbStatistics::BivariatePoisson::rv(theta1->getValue(), theta2->getValue(), theta0->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void BivariatePoissonDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == theta1)
    {
        theta1 = static_cast<const TypedDagNode<double>* >( newP );
    } else if (oldP == theta2)
    {
        theta2 = static_cast<const TypedDagNode<double>* >( newP );
    } else if (oldP == theta0)
    {
      theta0 = static_cast<const TypedDagNode<double>* >( newP );
    }
}
