#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_gamma.h"
#include "GammaDistribution.h"
#include "Real.h"
#include "RealPos.h"

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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_gamma::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_gamma::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Gamma probability distribution for positive real numbers.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_gamma::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The gamma distribution is the probability of the sum of exponentially distributed variables.";
    details += " Thus, it provides a natural prior distribution for parameters that could be considered as sums of exponential variables.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_gamma::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# lets simulate\n";
    example += "a <- rgamma(1000,shape=4,rate=4)\n";
    example += "# we expect a mean of 1\n";
    example += "mean(a)\n";
    example += "\n";
    example += "# create a random variable\n";
    example += "x ~ dnGamma(shape=4,rate=1)\n";
    example += "x\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_gamma::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_gamma::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnExponential" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_gamma::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Gamma Distribution";
    
    return title;
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
