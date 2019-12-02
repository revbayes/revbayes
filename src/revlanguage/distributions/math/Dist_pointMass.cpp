#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_pointMass.h"
#include "Real.h"
#include "RealPos.h"
#include "PointMassDistribution.h"

using namespace RevLanguage;

Dist_pointMass::Dist_pointMass() : ContinuousDistribution()
{
    
}


Dist_pointMass::~Dist_pointMass() {
    
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



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_pointMass::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Walker Pett" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_pointMass::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Point mass distribution.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_pointMass::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The point mass distribution, or Dirac delta function, has density f(x) = 1 when x is equal to the point mass value.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_pointMass::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "u ~ dnPointMass(1.2)\n";
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_pointMass::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_pointMass::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    // see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_pointMass::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Point Mass Distribution";
    
    return title;
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
