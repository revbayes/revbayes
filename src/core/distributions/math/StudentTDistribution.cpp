

#include "StudentTDistribution.h"

#include "DistributionStudentT.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** Student's T distribution Constructor
 * @param df Degrees of freedom
 *
 */

StudentTDistribution::StudentTDistribution(const TypedDagNode<std::int64_t> *df) : ContinuousDistribution( new double( 0.0 ) ),
degrees( df )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( degrees );
    
    *value = RbStatistics::StudentT::rv(degrees->getValue(), *GLOBAL_RNG);
}


StudentTDistribution::~StudentTDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double StudentTDistribution::cdf( void ) const
{
    return RbStatistics::StudentT::cdf(degrees->getValue(), *value );
}


StudentTDistribution* StudentTDistribution::clone( void ) const
{
    return new StudentTDistribution( *this );
}


double StudentTDistribution::computeLnProbability( void )
{
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::StudentT::lnPdf(degrees->getValue(), v);
}


double StudentTDistribution::getMax( void ) const
{
    return RbConstants::Double::inf;
}


double StudentTDistribution::getMin( void ) const
{
    return 0.0;
}


double StudentTDistribution::quantile(double p) const
{
    return RbStatistics::StudentT::quantile(degrees->getValue(), p);
}


void StudentTDistribution::redrawValue( void )
{
    *value = RbStatistics::StudentT::rv(degrees->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void StudentTDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == degrees)
    {
        degrees = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}


