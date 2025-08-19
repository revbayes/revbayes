#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_whiteNoise.h"
#include "RealPos.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlPositiveContinuousDistribution.h"
#include "TypeSpec.h"
#include "WhiteNoiseDistribution.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_whiteNoise::Dist_whiteNoise() : PositiveContinuousDistribution()
{
    
}


Dist_whiteNoise::~Dist_whiteNoise()
{
    
}



Dist_whiteNoise* Dist_whiteNoise::clone( void ) const
{
    
    return new Dist_whiteNoise(*this);
}


RevBayesCore::WhiteNoiseDistribution* Dist_whiteNoise::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* m   = static_cast<const RealPos &>( mu->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s   = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* t   = static_cast<const RealPos &>( time->getRevObject() ).getDagNode();
    RevBayesCore::WhiteNoiseDistribution* d = new RevBayesCore::WhiteNoiseDistribution(m, s, t);
    
    return d;
}


/* Get Rev type of object */
const std::string& Dist_whiteNoise::getClassType(void)
{
    
    static std::string rev_type = "Dist_whiteNoise";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Dist_whiteNoise::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( PositiveContinuousDistribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_whiteNoise::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "WhiteNoise";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_whiteNoise::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "mu",    RealPos::getClassTypeSpec(), "The mean of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "sigma", RealPos::getClassTypeSpec(), "The standard deviation of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "time" , RealPos::getClassTypeSpec(), "The time that the process has run.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_whiteNoise::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_whiteNoise::printValue(std::ostream& o) const
{
    
    o << "WhiteNoise(mu=";
    if ( mu != NULL )
    {
        o << mu->getName();
    }
    else
    {
        o << "?";
    }
    o << ", sigma=";
    if ( sigma != NULL )
    {
        o << sigma->getName();
    }
    else
    {
        o << "?";
    }
    o << ", time=";
    if ( time != NULL )
    {
        o << time->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_whiteNoise::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "time" )
    {
        time = var;
    }
    else
    {
        PositiveContinuousDistribution::setConstParameter(name, var);
    }
}
