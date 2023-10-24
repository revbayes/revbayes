#include "LinearDemographicFunction.h"

#include <cmath>

#include "Cloneable.h"
#include "RbException.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/**
 * @param[in]    N0   Pointer to the population size at the beginning of the linear period (towards the present)
 * @param[in]    N1   Pointer to the population size at the end of the linear period (towards the past)
 * @param[in]    t0   Pointer to the time at which the linear process started
 * @param[in]    t1   Pointer to the time of the beginning of the linear period (towards the present)
 */
LinearDemographicFunction::LinearDemographicFunction(const TypedDagNode<double>* N0, const TypedDagNode<double>* N1, const TypedDagNode<double>* t0, const TypedDagNode<double>* t1) : DemographicFunction(),
    theta_ancient( N1 ),
    theta_recent( N0 ),
    time_ancient( t1 ),
    time_recent( t0 )
{
    addVariable( N0 );
    addVariable( N1 );
    addVariable( t0 );
    addVariable( t1 );
}

/**
 * @param[in]    f    The linear demographic function to copy.
 */
LinearDemographicFunction::LinearDemographicFunction(const LinearDemographicFunction &f) : DemographicFunction(f),
    theta_ancient( f.theta_ancient ),
    theta_recent( f.theta_recent ),
    time_ancient( f.time_ancient ),
    time_recent( f.time_recent )
{
    
}


LinearDemographicFunction::~LinearDemographicFunction( void )
{
    
}

/**
 * @param[in]    f    The linear demographic function to copy.
 */
LinearDemographicFunction& LinearDemographicFunction::operator=(const LinearDemographicFunction &f)
{
    DemographicFunction::operator=( f );
    
    if ( this != &f )
    {
        theta_ancient   = f.theta_ancient;
        theta_recent    = f.theta_recent;
        time_ancient    = f.time_ancient;
        time_recent     = f.time_recent;
    }
    
    return *this;
}

/**
 * This is similar to the copy constructor but useful in inheritance.
 *
 * The clone function is a convenience function to create proper copies of inherited objects.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * @return A new copy of myself
 */
LinearDemographicFunction* LinearDemographicFunction::clone( void ) const
{
    
    return new LinearDemographicFunction(*this);
}

/**
 * @param[in]   t    Time
 *
 * @return  N(t)
 */
double LinearDemographicFunction::getDemographic(double t) const
{
    double N0 = theta_recent->getValue();
    double N1 = theta_ancient->getValue();
    double t0 = time_recent->getValue();
    double t1 = time_ancient->getValue();
    
    if ( t1 < t0 || t0 < 0 || N1 < 0 || t < t0 || t > t1)
    {
        throw RbException("Impossible parameter values in Linear growth/decline demographic functions.");
    }
    
    if ( N0 == N1 )
    {
        return N0;
    }
    else
    {
        double alpha = ( N1-N0 ) / (t1 - t0);
        return N0 + (t-t0) * alpha;
    }
    
}

/**
 * @param[in]   start   Time at which the interval starts
 * @param[in]   finish  Time at which the interval ends
 *
 * @return  Integral 1/N(x) dx between start and finish.
 */
double LinearDemographicFunction::getIntegral(double start, double finish) const
{
    double N0 = theta_recent->getValue();
    double N1 = theta_ancient->getValue();
    double t0 = time_recent->getValue();
    double t1 = time_ancient->getValue();
    
    if ( t1 < t0 || t0 < 0 || N1 < 0 || start < t0 || start > t1 || finish < t0 || finish > t1 )
    {
        throw RbException("Impossible parameter values in Linear growth/decline demographic functions.");
    }
    
    
    if ( N0 == N1 )
    {
        double delta = finish - start;
        return delta / N0;
    }
    else
    {
        double alpha = ( N1-N0 ) / (t1 - t0);
        return ( log( N0 + (finish-t0) * alpha ) - log( N0 + (start-t0) * alpha ) ) / alpha;
    }
    
}

/**
 * @param[in]   time    Current time in coalescent simulation process
 * @param[in]   lambda
 *
 * @return Waiting Time until next coalescent event
 */
double LinearDemographicFunction::getWaitingTime(double time, double lambda) const
{
    double N0 = theta_recent->getValue();
    double N1 = theta_ancient->getValue();
    double t0 = time_recent->getValue();
    double t1 = time_ancient->getValue();
    
    if ( t1 < t0 || t0 < 0 || N1 < 0 || time < t0 || time > t1)
    {
        throw RbException("Impossible parameter values in Linear growth/decline demographic functions.");
    }
    
    if ( N0 == N1 )
    {
        return N0 * lambda;
    }
    else
    {
        double alpha = ( N1-N0 ) / (t1 - t0);
        return (N0 + (time-t0) * alpha) * lambda;
    }
}

/**
 * @param[in]   old_node    Pointer to the DAG node to be replaced
 * @param[in]   new_node    Pointer to the DAG node replacing the other
 */
void LinearDemographicFunction::swapNodeInternal(const DagNode *old_node, const DagNode *new_node)
{
    
    if (old_node == theta_ancient)
    {
        theta_ancient = static_cast<const TypedDagNode<double>* >( new_node );
    }
    
    if (old_node == theta_recent)
    {
        theta_recent = static_cast<const TypedDagNode<double>* >( new_node );
    }
    
    if (old_node == time_ancient)
    {
        time_ancient = static_cast<const TypedDagNode<double>* >( new_node );
    }
    
    if (old_node == time_recent)
    {
        time_recent = static_cast<const TypedDagNode<double>* >( new_node );
    }
    
}


std::ostream& operator<<(std::ostream& o, const LinearDemographicFunction& x)
{
    return o;
}
