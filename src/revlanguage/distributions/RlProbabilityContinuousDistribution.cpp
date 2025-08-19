#include "RlProbabilityContinuousDistribution.h"

#include "RlContinuousStochasticNode.h"
#include "TypedDistribution.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

using namespace RevLanguage;

ProbabilityContinuousDistribution::ProbabilityContinuousDistribution() : TypedDistribution<Probability>()
{
    
}



ProbabilityContinuousDistribution::ProbabilityContinuousDistribution( const ProbabilityContinuousDistribution &d ) : TypedDistribution<Probability>(d)
{
    
}



ProbabilityContinuousDistribution::~ProbabilityContinuousDistribution()
{
    
}



Probability* ProbabilityContinuousDistribution::createRandomVariable(void) const
{
    
    RevBayesCore::ContinuousDistribution* d = createDistribution();
    RevBayesCore::TypedDagNode<double>* rv  = new ContinuousStochasticNode("", d, this->clone() );
    
    return new Probability(rv);
}



/* Get Rev type of object */
const std::string& ProbabilityContinuousDistribution::getClassType(void)
{
    
    static std::string rev_type = "ProbabilityContinuousDistribution";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& ProbabilityContinuousDistribution::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Probability>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}
