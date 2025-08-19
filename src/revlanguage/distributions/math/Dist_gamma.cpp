#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_gamma.h"
#include "GammaDistribution.h"
#include "RealPos.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlPositiveContinuousDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_gamma::Dist_gamma() : PositiveContinuousDistribution()
{
    
}


Dist_gamma::~Dist_gamma()
{
    
}



Dist_gamma* Dist_gamma::clone( void ) const
{
    
    return new Dist_gamma(*this);
}


RevBayesCore::GammaDistribution* Dist_gamma::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* sh  = static_cast<const RealPos &>( shape->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* r   = static_cast<const RealPos &>( rate->getRevObject() ).getDagNode();
    RevBayesCore::GammaDistribution* d      = new RevBayesCore::GammaDistribution(sh, r);
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_gamma::getClassType(void)
{
    
    static std::string rev_type = "Dist_gamma";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_gamma::getClassTypeSpec(void)
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
std::string Dist_gamma::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "gamma";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_gamma::getParameterRules(void) const
{
    
    static MemberRules distGammaMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        distGammaMemberRules.push_back( new ArgumentRule( "shape", RealPos::getClassTypeSpec(), "The shape parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distGammaMemberRules.push_back( new ArgumentRule( "rate" , RealPos::getClassTypeSpec(), "The rate parameter (rate = 1/scale).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distGammaMemberRules;
}


const TypeSpec& Dist_gamma::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_gamma::printValue(std::ostream& o) const {
    
    o << "gamma(shape=";
    if ( shape != NULL ) 
    {
        o << shape->getName();
    } 
    else 
    {
        o << "?";
    }
    o << ", rate=";
    if ( rate != NULL ) 
    {
        o << rate->getName();
    } 
    else 
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_gamma::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "shape" ) 
    {
        shape = var;
    }
    else if ( name == "rate" ) 
    {
        rate = var;
    }
    else  
    {
        PositiveContinuousDistribution::setConstParameter(name, var);
    }
}
