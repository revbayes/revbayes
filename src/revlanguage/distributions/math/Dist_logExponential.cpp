#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_logExponential.h"
#include "LogExponentialDistribution.h"
#include "Real.h"
#include "RealPos.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "RlDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_logExponential::Dist_logExponential() : ContinuousDistribution()
{
    
}


Dist_logExponential::~Dist_logExponential()
{
    
}



Dist_logExponential* Dist_logExponential::clone( void ) const
{
    return new Dist_logExponential(*this);
}


RevBayesCore::LogExponentialDistribution* Dist_logExponential::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>*       l = static_cast<const Real &>( lambda->getRevObject() ).getDagNode();
    RevBayesCore::LogExponentialDistribution* d = new RevBayesCore::LogExponentialDistribution( l );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_logExponential::getClassType(void)
{
    
    static std::string rev_type = "Dist_logExponential";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_logExponential::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_logExponential::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "logExponential";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_logExponential::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The rate of the exonential distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_logExponential::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_logExponential::printValue(std::ostream& o) const
{
    
    o << " logExponential(min=";
    if ( lambda != NULL )
    {
        o << lambda->getName();
    }
    else
    {
        o << "?";
    }
    
    o << ")";
}


/** Set a member variable */
void Dist_logExponential::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "lambda" )
    {
        lambda = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}

