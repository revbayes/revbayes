#include "OrnsteinUhlenbeckProcess.h"

#include <cmath>

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
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
 * @param   x    The starting value for the process
 * @param   t    The equilibrium value of the process
 * @param   a    The strength of the spring
 * @param   s    The standard deviation
 * @param   ti   The time during which the process is applied
 */
OrnsteinUhlenbeckProcess::OrnsteinUhlenbeckProcess(const TypedDagNode<double> *x, const TypedDagNode<double> *t, const TypedDagNode<double> *a, const TypedDagNode<double> *s, const TypedDagNode<double> *ti) : ContinuousDistribution( new double( 0.0 ) ),
    x0( x ),
    theta( t ),
    alpha( a ),
    sigma( s ),
    time( ti )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( x0 );
    addParameter( theta );
    addParameter( alpha );
    addParameter( sigma );
    addParameter( time );
    
    redrawValue();
}


/**
 * Get the cumulative density of the current value.
 *
 * \return    The cumulative density.
 */
double OrnsteinUhlenbeckProcess::cdf( void ) const
{
    double a = alpha->getValue();
    double t = time->getValue();
    double e = exp(-a*t);
    double e2 = exp(-2.0*a*t);
    double m = e * x0->getValue() + (1.0-e)*theta->getValue();
    double s = sigma->getValue() * sqrt( (1-e2) / 2 / a );
    
    return RbStatistics::Normal::cdf( m, s, *value);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
OrnsteinUhlenbeckProcess* OrnsteinUhlenbeckProcess::clone( void ) const
{
    
    return new OrnsteinUhlenbeckProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return   The log-transformed probability density.
 */
double OrnsteinUhlenbeckProcess::computeLnProbability( void )
{
    double a = alpha->getValue();
    double t = time->getValue();
    double e = exp(-a*t);
    double e2 = exp(-2.0*a*t);
    double m = e * x0->getValue() + (1.0-e)*theta->getValue();
    double s = sigma->getValue() * sqrt( (1-e2) / 2 / a );
    
    return RbStatistics::Normal::lnPdf( m, s, *value);
}


/**
 * Get the maximum a value drawn from this distribution can take.
 *
 * \return    Positive infinity as the maximum value.
 */
double OrnsteinUhlenbeckProcess::getMax( void ) const
{
    return RbConstants::Double::inf;
}


/**
 * Get the mininum a value drawn from this distribution can take.
 *
 * \return    Negative infinity as the minimum value.
 */
double OrnsteinUhlenbeckProcess::getMin( void ) const
{
    return RbConstants::Double::neginf;
}


/**
 * Get the quantile for the probability p.
 *
 * \param[in]   p   The probability for which the quantile should be computed.
 *
 * \return    The quantile.
 *
 * @todo This is not yet implemented
 */
double OrnsteinUhlenbeckProcess::quantile(double p) const
{
    throw RbException("Quantile function of Ornstein-Uhlenbeck process not implemented.");
}


/**
 * Redrawing a new value from the process and storing it as a member.
 */
void OrnsteinUhlenbeckProcess::redrawValue( void )
{
    double a = alpha->getValue();
    double t = time->getValue();
    double e = exp(-a*t);
    double e2 = exp(-2.0*a*t);
    double m = e * x0->getValue() + (1.0-e)*theta->getValue();
    double s = sigma->getValue() * sqrt( (1-e2) / 2 / a );
    
    *value = RbStatistics::Normal::rv( m, s, *GLOBAL_RNG);
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void OrnsteinUhlenbeckProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == x0)
    {
        x0 = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == theta)
    {
        theta = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == sigma)
    {
        sigma = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == time)
    {
        time = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}
