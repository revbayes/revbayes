#include "RlPositiveContinuousDistribution.h"

#include "RlContinuousStochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

PositiveContinuousDistribution::PositiveContinuousDistribution() : TypedDistribution<RealPos>()
{
    
}



PositiveContinuousDistribution::PositiveContinuousDistribution( const PositiveContinuousDistribution &d ) : TypedDistribution<RealPos>(d)
{
    
}



PositiveContinuousDistribution::~PositiveContinuousDistribution()
{
    
}



RealPos* PositiveContinuousDistribution::createRandomVariable(void) const
{
    
    RevBayesCore::ContinuousDistribution* d = createDistribution();
    RevBayesCore::TypedDagNode<double>* rv  = new ContinuousStochasticNode("", d, this->clone() );
    
    return new RealPos(rv);
}



/* Get Rev type of object */
const std::string& PositiveContinuousDistribution::getClassType(void)
{
    
    static std::string rev_type = "PositiveContinuousDistribution";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& PositiveContinuousDistribution::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<RealPos>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}

