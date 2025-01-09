#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_UniformInteger.h"
#include "UniformIntegerDistribution.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "Integer.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;

Dist_UniformInteger::Dist_UniformInteger() : TypedDistribution<Integer>()
{
    
}


Dist_UniformInteger::~Dist_UniformInteger()
{
    
}



Dist_UniformInteger* Dist_UniformInteger::clone( void ) const
{
    return new Dist_UniformInteger(*this);
}


RevBayesCore::UniformIntegerDistribution* Dist_UniformInteger::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<long>* l   = static_cast<const Integer &>( lower->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<long>* u   = static_cast<const Integer &>( upper->getRevObject() ).getDagNode();
    RevBayesCore::UniformIntegerDistribution* d    = new RevBayesCore::UniformIntegerDistribution(l, u);
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_UniformInteger::getClassType(void)
{
    
    static std::string rev_type = "Dist_UniformInteger";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_UniformInteger::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Integer>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_UniformInteger::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "unifInt" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_UniformInteger::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "UniformInteger";
    
    return d_name;
}


/** Return member rules */
const MemberRules& Dist_UniformInteger::getParameterRules(void) const
{
    
    static MemberRules distUnifMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        distUnifMemberRules.push_back( new ArgumentRule( "lower", Integer::getClassTypeSpec(), "The lower bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distUnifMemberRules.push_back( new ArgumentRule( "upper", Integer::getClassTypeSpec(), "The upper bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distUnifMemberRules;
}


const TypeSpec& Dist_UniformInteger::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_UniformInteger::printValue(std::ostream& o) const
{
    
    o << " unif (lower=";
    if ( lower != NULL )
    {
        o << lower->getName();
    }
    else
    {
        o << "?";
    }
    o << ", upper=";
    if ( lower != NULL )
    {
        o << lower->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_UniformInteger::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "lower" )
    {
        lower = var;
    }
    else if ( name == "upper" )
    {
        upper = var;
    }
    else
    {
        TypedDistribution<Integer>::setConstParameter(name, var);
    }
}
