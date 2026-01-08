#include "SoftBoundUniformNormalDistribution.h"

#include <cmath>

#include "DistributionNormal.h"
#include "DistributionUniform.h"
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
 * \param[in]   mi    The min of the uniform distribution.
 * \param[in]   ma    The max of the uniform distribution.
 * \param[in]   sd    The standard deviation of the normal distribution.
 * \param[in]   p     The probability that the realization came from the uniform distribution.
 */
SoftBoundUniformNormalDistribution::SoftBoundUniformNormalDistribution(const TypedDagNode<double> *mi, const TypedDagNode<double> *ma, const TypedDagNode<double> *sd, const TypedDagNode<double> *p, BOUNDS b) : ContinuousDistribution( new double( 0.0 ) ),
    min( mi ),
    max( ma ),
    stDev( sd ),
    prob( p ),
    soft_bound( b )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( min );
    addParameter( max );
    addParameter( stDev );
    addParameter( prob );
    
    redrawValue();
}


/**
 * Get the cumulative density of the current value.
 *
 * \return    The cumulative density.
 */
double SoftBoundUniformNormalDistribution::cdf( void ) const
{

    throw RbException("Missing implementation.");
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
SoftBoundUniformNormalDistribution* SoftBoundUniformNormalDistribution::clone( void ) const
{
    
    return new SoftBoundUniformNormalDistribution( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return   The log-transformed probability density.
 */
double SoftBoundUniformNormalDistribution::computeLnProbability( void )
{
    double side_factor = ( soft_bound == BOTH ? 2.0 : 1.0 );       // This is used for the computation of the sd or p if only one is provided.

    double min_val  = min->getValue();
    double max_val  = max->getValue();
    double sd       = 1.0;
    double p        = 0.95;
    
    if ( stDev == NULL && prob == NULL )
    {
        throw RbException( "Cannot compute sd and prob in SoftBoundUniformNormal distribution if neither is provided." );
    }
    else if ( stDev == NULL )
    {
        p   = prob->getValue();
        sd  = (1-p) * (max_val-min_val) / (side_factor*p*RbConstants::SQRT_2PI);
    }
    else if ( prob == NULL )
    {
        sd  = stDev->getValue();
        p   = 1 / ( side_factor/(max_val-min_val) * RbConstants::SQRT_2PI * sd + 1 );
    }
    else
    {
        p   = prob->getValue();
        sd  = stDev->getValue();
    }
    
    if ( *value < min_val )
    {
        if ( soft_bound == BOTH )
        {
            return log((1-p)/2.0) + RbStatistics::Normal::lnPdf( 0.0, sd, *value - min_val );
        }
        else if ( soft_bound == LOWER )
        {
            return log(1-p) + RbStatistics::Normal::lnPdf( 0.0, sd, *value - min_val );
        }
        else
        {
            return RbConstants::Double::neginf;
        }
    }
    else if ( *value > max_val )
    {
        
        if ( soft_bound == BOTH )
        {
            return log((1-p)/2.0) + RbStatistics::Normal::lnPdf( 0.0, sd, *value - max_val );
        }
        else if ( soft_bound == UPPER )
        {
            return log(1-p) + RbStatistics::Normal::lnPdf( 0.0, sd, *value - max_val );
        }
        else
        {
            return RbConstants::Double::neginf;
        }

    }
    else
    {
        return log( p ) + RbStatistics::Uniform::lnPdf( min_val, max_val, *value);
    }

}


/**
 * Get the maximum a value drawn from this distribution can take.
 *
 * \return    Positive infinity as the maximum value.
 */
double SoftBoundUniformNormalDistribution::getMax( void ) const
{
    return RbConstants::Double::inf;
}


/**
 * Get the mininum a value drawn from this distribution can take.
 *
 * \return    Negative infinity as the minimum value.
 */
double SoftBoundUniformNormalDistribution::getMin( void ) const
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
double SoftBoundUniformNormalDistribution::quantile(double p) const
{
    throw RbException("Quantile function of SoftBoundUniform normal distribution not implemented.");
}


/**
 * Redrawing a new value from the process and storing it as a member.
 */
void SoftBoundUniformNormalDistribution::redrawValue( void )
{
    double side_factor = ( soft_bound == BOTH ? 2.0 : 1.0 );       // This is used for the computation of the sd or p if only one is provided.

    double min_val  = min->getValue();
    double max_val  = max->getValue();
    double sd       = 1.0;
    double p        = 0.95;
    
    if ( stDev == NULL && prob == NULL )
    {
        throw RbException( "Cannot compute sd and prob in SoftBoundUniformNormal distribution if neither is provided." );
    }
    else if ( stDev == NULL )
    {
        p   = prob->getValue();
        sd  = (1-p) * (max_val-min_val) / (side_factor*p*RbConstants::SQRT_2PI);
    }
    else if ( prob == NULL )
    {
        sd  = stDev->getValue();
        p   = 1 / ( side_factor/(max_val-min_val) * RbConstants::SQRT_2PI * sd + 1 );
    }
    else
    {
        p   = prob->getValue();
        sd  = stDev->getValue();
    }
    
    double u = GLOBAL_RNG->uniform01();
    if ( u < p )
    {
        *value = RbStatistics::Uniform::rv(min_val, max_val, *GLOBAL_RNG);
    }
    else
    {
        double x = RbStatistics::Normal::rv(0.0, sd, *GLOBAL_RNG);
        if ( soft_bound == LOWER )
        {
            x = fmin(x,-x);
        }
        else if ( soft_bound == UPPER )
        {
            x = fmax(x,-x);
        }
        if ( x > 0.0 )
        {
            *value = max_val + x;
        }
        else
        {
            *value = min_val + x;
        }
        
    }
    
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void SoftBoundUniformNormalDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == min)
    {
        min = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == max)
    {
        max = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == stDev)
    {
        stDev = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == prob)
    {
        prob = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}
