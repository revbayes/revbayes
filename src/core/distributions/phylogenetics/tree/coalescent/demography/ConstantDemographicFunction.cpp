#include "ConstantDemographicFunction.h"

#include "TypedDagNode.h"
#include "Cloneable.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


ConstantDemographicFunction::ConstantDemographicFunction(const TypedDagNode<double>* t) : DemographicFunction(),
    theta( t )
{
    // add the parameter
    addVariable( theta );
    
}


ConstantDemographicFunction::ConstantDemographicFunction(const ConstantDemographicFunction &f) : DemographicFunction(f),
    theta( f.theta )
{
    
}


ConstantDemographicFunction::~ConstantDemographicFunction( void )
{
    
}


ConstantDemographicFunction& ConstantDemographicFunction::operator=(const ConstantDemographicFunction &f)
{
    DemographicFunction::operator=( f );
    
    if ( this != &f )
    {
        theta = f.theta;
    }
    
    return *this;
}


ConstantDemographicFunction* ConstantDemographicFunction::clone( void ) const
{
    
    return new ConstantDemographicFunction(*this);
}


double ConstantDemographicFunction::getDemographic(double t) const
{
    
    return theta->getValue();
}


double ConstantDemographicFunction::getIntegral(double start, double finish) const
{
    double delta = finish - start;
    return delta / theta->getValue();
}

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
