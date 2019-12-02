#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_logExponential.h"
#include "LogExponentialDistribution.h"
#include "Real.h"
#include "RealPos.h"

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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_logExponential::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_logExponential::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A real number x has a log-Exponential distribution if y = exp(x) has Exponential distribution.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_logExponential::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The log-Exponential distribution is defined over real numbers.";
    details += " Saying that x is log-Exponential is equivalent to saying that y = exp(x) is Exponential.";
    details += " The log-Exponential distribution therefore expresses lack of information about the order of magnitude of a scale parameter: ";
    details += " if x has a log-Exponential distribution, then it has equal chance to be contained by any of the intervals of the form (10^k, 10^(k+1)) within the allowed range.\n\n";
    
    details += "The density is p(x) = ???, which can be seen by defining x = ln(y) where y has Exponential distribution and apply the change-of-variable formula.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_logExponential::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# a log-Exponential prior over the rate of change of a Brownian trait (or a Brownian relaxed clock)\n";
    example += "trueTree = readTrees(\"data/primates.tree\")[1]\n";
    example += "log_sigma ~ dnLogExponential(lambda=1)\n";
    example += "sigma := exp(log_sigma)\n";
    example += "X ~ dnBrownian(trueTree,sigma)\n";
    example += "# ...\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_logExponential::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_logExponential::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnExponential" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_logExponential::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Log-Exponential Distribution";
    
    return title;
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

