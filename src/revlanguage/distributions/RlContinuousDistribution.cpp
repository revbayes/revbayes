#include "RlContinuousDistribution.h"

#include "RlContinuousStochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

ContinuousDistribution::ContinuousDistribution() : TypedDistribution<Real>()
{
    
}



ContinuousDistribution::ContinuousDistribution( const ContinuousDistribution &d ) : TypedDistribution<Real>(d)
{
    
}



ContinuousDistribution::~ContinuousDistribution()
{
    
}



Real* ContinuousDistribution::createRandomVariable(void) const
{
    
    RevBayesCore::ContinuousDistribution* d = createDistribution();
    RevBayesCore::TypedDagNode<double>* rv  = new ContinuousStochasticNode("", d, this->clone() );
    
    return new Real(rv);
}



/* Get Rev type of object */
const std::string& ContinuousDistribution::getClassType(void)
{
    
    static std::string rev_type = "ContinuousDistribution";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& ContinuousDistribution::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Real>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}

