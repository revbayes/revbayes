#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_pointMass.h"
#include "Real.h"
#include "PointMassDistribution.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_pointMass::Dist_pointMass() : ContinuousDistribution()
{
    
}


Dist_pointMass::~Dist_pointMass()
{
    
}



Dist_pointMass* Dist_pointMass::clone( void ) const
{
    return new Dist_pointMass (*this);
}


RevBayesCore::PointMassDistribution* Dist_pointMass::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* l   = static_cast<const Real &>( value->getRevObject() ).getDagNode();
    RevBayesCore::PointMassDistribution* d    = new RevBayesCore::PointMassDistribution(l);
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_pointMass::getClassType(void)
{
    
    static std::string rev_type = "Dist_pointMass";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_pointMass::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_pointMass::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "DiracDelta" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_pointMass::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PointMass";
    
    return d_name;
}


/** Return member rules */
const MemberRules& Dist_pointMass::getParameterRules(void) const
{
    
    static MemberRules distUnifMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        distUnifMemberRules.push_back( new ArgumentRule( "value", Real::getClassTypeSpec(), "The value where mass is concentrated.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distUnifMemberRules;
}


const TypeSpec& Dist_pointMass::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_pointMass::printValue(std::ostream& o) const
{
    
    o << " pointMass (value=";
    if ( value != NULL ) {
        o << value->getName();
    } else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_pointMass::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "value" )
    {
        value = var;
    }
    else 
    {
        ContinuousDistribution::setConstParameter(name, var);
    }
}
