#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_lnormNegativeOffset.h"
#include "LognormalDistribution.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_lnormNegativeOffset::Dist_lnormNegativeOffset() : ContinuousDistribution()
{
    
}



Dist_lnormNegativeOffset* Dist_lnormNegativeOffset::clone( void ) const
{
    return new Dist_lnormNegativeOffset(*this);
}


RevBayesCore::LognormalDistribution* Dist_lnormNegativeOffset::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* m   = static_cast<const Real &>( mean->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s   = static_cast<const RealPos &>( sd->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* o   = static_cast<const Real &>( offset->getRevObject() ).getDagNode();
    RevBayesCore::LognormalDistribution* d  = new RevBayesCore::LognormalDistribution(m, s, o);
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_lnormNegativeOffset::getClassType(void) {
    
    static std::string rev_type = "Dist_lnormNegativeOffset";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_lnormNegativeOffset::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_lnormNegativeOffset::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "lnorm" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_lnormNegativeOffset::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "lognormal";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_lnormNegativeOffset::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        memberRules.push_back( new ArgumentRule( "mean",   Real::getClassTypeSpec()   , "The mean in log-space (observed mean is exp(m))." , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "sd"  ,   RealPos::getClassTypeSpec(), "The standard deviation in log-space."             , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "offset", Real::getClassTypeSpec()   , "The offset in real-space."                        , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        
        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_lnormNegativeOffset::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_lnormNegativeOffset::printValue(std::ostream& o) const {
    
    o << "lognormal(mean=";
    if ( mean != NULL ) {
        o << mean->getName();
    } else {
        o << "?";
    }
    o << ", sd=";
    if ( sd != NULL ) {
        o << sd->getName();
    } else {
        o << "?";
    }
    o << ", offset=";
    if ( offset != NULL ) {
        o << offset->getName();
    } else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_lnormNegativeOffset::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "mean" ) 
    {
        mean = var;
    }
    else if ( name == "sd" ) 
    {
        sd = var;
    }
    else if ( name == "offset" )
    {
        offset = var;
    }
    else 
    {
        ContinuousDistribution::setConstParameter(name, var);
    }
}
