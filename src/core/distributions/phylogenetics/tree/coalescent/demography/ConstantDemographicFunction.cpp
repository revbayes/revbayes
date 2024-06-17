#include "ConstantDemographicFunction.h"

#include "TypedDagNode.h"
#include "Cloneable.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

/**
 * @param[in]   t    Pointer to theta, the effective population size. Here: theta is constant.
 */
ConstantDemographicFunction::ConstantDemographicFunction(const TypedDagNode<double>* t) : DemographicFunction(),
    theta( t )
{
    // add the parameter
    addVariable( theta );
    
}

/**
 * @param[in]    f    The constant demographic function to copy.
 */
ConstantDemographicFunction::ConstantDemographicFunction(const ConstantDemographicFunction &f) : DemographicFunction(f),
    theta( f.theta )
{
    
}


ConstantDemographicFunction::~ConstantDemographicFunction( void )
{
    
}


/**
 * @param[in]    f    The constant demographic function to copy.
 */
ConstantDemographicFunction& ConstantDemographicFunction::operator=(const ConstantDemographicFunction &f)
{
    DemographicFunction::operator=( f );
    
    if ( this != &f )
    {
        theta = f.theta;
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
ConstantDemographicFunction* ConstantDemographicFunction::clone( void ) const
{
    
    return new ConstantDemographicFunction(*this);
}

/**
 * @param[in]   t    Time
 *
 * @return  N(t)
 */
double ConstantDemographicFunction::getDemographic(double t) const
{
    
    return theta->getValue();
}

/**
 * @param[in]   start   Time at which the interval starts
 * @param[in]   finish  Time at which the interval ends
 *
 * @return  Integral 1/N(x) dx between start and finish.
 */
double ConstantDemographicFunction::getIntegral(double start, double finish) const
{
    double delta = finish - start;
    return delta / theta->getValue();
}

/**
 * @param[in]   time    Current time in coalescent simulation process
 * @param[in]   lambda  Waiting time under a standardized coalescent with constant population of theta=1
 *
 * @return Waiting Time until next coalescent event
 */
double ConstantDemographicFunction::getWaitingTime(double time, double lambda, double ploidy_factor) const
{
    return theta->getValue() * lambda * ploidy_factor;
}

/**
 * @param[in]   old_node    Pointer to the DAG node to be replaced
 * @param[in]   new_node    Pointer to the DAG node replacing the other
 */
void ConstantDemographicFunction::swapNodeInternal(const DagNode *old_node, const DagNode *new_node)
{

    if (old_node == theta)
    {
        theta = static_cast<const TypedDagNode<double>* >( new_node );
    }
    
}


std::ostream& operator<<(std::ostream& o, const ConstantDemographicFunction& x)
{
    return o;
}
