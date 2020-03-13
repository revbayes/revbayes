#include "BimodalLognormalDistribution.h"

#include <cmath>

#include "DistributionLognormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Default constructor.
 * The default constructor does nothing except allocating the object.
 *
 * \param[in]   m1    The mean of the first lognormal distribution.
 * \param[in]   m2    The mean of the second lognormal distribution.
 * \param[in]   s1    The standard deviation of the first lognormal distribution.
 * \param[in]   s2    The standard deviation of the second lognormal distribution.
 * \param[in]   p     The probability that the realization came from the first lognormal distribution.
 */
BimodalLognormalDistribution::BimodalLognormalDistribution(const TypedDagNode<double> *m1, const TypedDagNode<double> *m2, const TypedDagNode<double> *s1, const TypedDagNode<double> *s2, const TypedDagNode<double> *p) : ContinuousDistribution( new double( 0.0 ) ),
    mean1( m1 ),
    mean2( m2 ),
    stDev1( s1 ),
    stDev2( s2 ),
    p( p )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( m1 );
    addParameter( m2 );
    addParameter( s1 );
    addParameter( s2 );
    addParameter( p );

    double u = GLOBAL_RNG->uniform01();
    if ( u < p->getValue() )
    {
        *value = RbStatistics::Lognormal::rv(mean1->getValue(), stDev1->getValue(), *GLOBAL_RNG);
    }
    else
    {
        *value = RbStatistics::Lognormal::rv(mean2->getValue(), stDev2->getValue(), *GLOBAL_RNG);
    }
}


/**
 * Get the cumulative density of the current value.
 *
 * \return    The cumulative density.
 */
double BimodalLognormalDistribution::cdf( void ) const
{
    return p->getValue() * RbStatistics::Lognormal::cdf( mean1->getValue(), stDev1->getValue(), *value) + (1.0 - p->getValue()) * RbStatistics::Lognormal::cdf( mean2->getValue(), stDev2->getValue(), *value);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
BimodalLognormalDistribution* BimodalLognormalDistribution::clone( void ) const
{

    return new BimodalLognormalDistribution( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return   The log-transformed probability density.
 */
double BimodalLognormalDistribution::computeLnProbability( void )
{
    return log(p->getValue() * RbStatistics::Lognormal::pdf( mean1->getValue(), stDev1->getValue(), *value) + (1.0 - p->getValue()) * RbStatistics::Lognormal::pdf( mean2->getValue(), stDev2->getValue(), *value) );
}


/**
 * Get the maximum a value drawn from this distribution can take.
 *
 * \return    Positive infinity as the maximum value.
 */
double BimodalLognormalDistribution::getMax( void ) const
{
    return RbConstants::Double::inf;
}


/**
 * Get the mininum a value drawn from this distribution can take.
 *
 * \return    Negative infinity as the minimum value.
 */
double BimodalLognormalDistribution::getMin( void ) const
{
    return RbConstants::Double::neginf;
}


/**
 * Get the quantile for the probability p.
 *
 * \param[in]   p   The probability for which the quantile should be computed.
 *
 * \return    The quantile.
 */
double BimodalLognormalDistribution::quantile(double p) const
{
    throw RbException("Quantile function of bimodal lognormal distribution not implemented.");
}


/**
 * Redrawing a new value from the process and storing it as a member.
 */
void BimodalLognormalDistribution::redrawValue( void )
{
    double u = GLOBAL_RNG->uniform01();
    if ( u < p->getValue() )
    {
        *value = RbStatistics::Lognormal::rv(mean1->getValue(), stDev1->getValue(), *GLOBAL_RNG);
    }
    else
    {
        *value = RbStatistics::Lognormal::rv(mean2->getValue(), stDev2->getValue(), *GLOBAL_RNG);
    }
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void BimodalLognormalDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == mean1)
    {
        mean1 = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == mean2)
    {
        mean2 = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == stDev1)
    {
        stDev1 = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == stDev2)
    {
        stDev2 = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == p)
    {
        p = static_cast<const TypedDagNode<double>* >( newP );
    }

}
