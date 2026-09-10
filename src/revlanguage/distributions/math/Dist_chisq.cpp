#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ChisqDistribution.h"
#include "Dist_chisq.h"
#include "Natural.h"
#include "RealPos.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDistribution.h"
#include "RlPositiveContinuousDistribution.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"
#include "TypedDistribution.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_chisq::Dist_chisq(void) : PositiveContinuousDistribution()
{
    
}


Dist_chisq::~Dist_chisq(void)
{
    
}



Dist_chisq* Dist_chisq::clone( void ) const
{
    
    return new Dist_chisq(*this);
}


RevBayesCore::ChisqDistribution* Dist_chisq::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<std::int64_t>*    df = static_cast<const Natural     &>( degrees->getRevObject() ).getDagNode();
    RevBayesCore::ChisqDistribution* d  = new RevBayesCore::ChisqDistribution( df );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_chisq::getClassType(void)
{
    
    static std::string rev_type = "Dist_chisq";
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_chisq::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_chisq::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "chisq";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_chisq::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "df", Natural::getClassTypeSpec(), "The degrees of freedom.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_chisq::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_chisq::printValue(std::ostream& o) const
{
    
    o << "chisq(df=";
    if ( degrees != NULL )
    {
        o << degrees->getName();
    }
    else
    {
        o << "?";
    }
    
    o << ")";
}


/** Set a member variable */
void Dist_chisq::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "df" )
    {
        degrees = var;
    }
    else
    {
        TypedDistribution<RealPos>::setConstParameter(name, var);
    }
}
