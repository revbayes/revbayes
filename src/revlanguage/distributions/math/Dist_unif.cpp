#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_unif.h"
#include "Real.h"
#include "RealPos.h"
#include "UniformDistribution.h"

using namespace RevLanguage;

Dist_unif::Dist_unif() : ContinuousDistribution()
{
    
}


Dist_unif::~Dist_unif() {
    
}



Dist_unif* Dist_unif::clone( void ) const
{
    return new Dist_unif (*this);
}


RevBayesCore::UniformDistribution* Dist_unif::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* l   = static_cast<const Real &>( lower->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* u   = static_cast<const Real &>( upper->getRevObject() ).getDagNode();
    RevBayesCore::UniformDistribution* d    = new RevBayesCore::UniformDistribution(l, u);
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_unif::getClassType(void)
{
    
    static std::string rev_type = "Dist_unif";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_unif::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_unif::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "unif" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_unif::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Uniform";
    
    return d_name;
}



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_unif::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael Landis" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_unif::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Uniform distribution with lower and uppper bounds.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_unif::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The uniform distribution has density, f(x) = 1/(b-a), where b is the upper bound and a is the lower bound.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_unif::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# set the lower bound\n";
    example += "a <- -2.5\n";
    example += "# set the upper bound\n";
    example += "b <- -3.9\n";
    example += "# create a stochastic node with a uniform prior\n";
    example += "u ~ dnUniform(a, b)\n";
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_unif::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_unif::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    // see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_unif::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Uniform Distribution";
    
    return title;
}


/** Return member rules */
const MemberRules& Dist_unif::getParameterRules(void) const
{
    
    static MemberRules distUnifMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        distUnifMemberRules.push_back( new ArgumentRule( "lower", Real::getClassTypeSpec(), "The lower bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distUnifMemberRules.push_back( new ArgumentRule( "upper", Real::getClassTypeSpec(), "The upper bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distUnifMemberRules;
}


const TypeSpec& Dist_unif::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_unif::printValue(std::ostream& o) const
{
    
    o << " unif (lower=";
    if ( lower != NULL ) {
        o << lower->getName();
    } else {
        o << "?";
    }
    o << ", upper=";
    if ( lower != NULL ) {
        o << lower->getName();
    } else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_unif::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
        ContinuousDistribution::setConstParameter(name, var);
    }
}
