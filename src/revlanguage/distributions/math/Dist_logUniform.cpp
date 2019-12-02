#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_logUniform.h"
#include "LogUniformDistribution.h"
#include "Real.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_logUniform::Dist_logUniform() : PositiveContinuousDistribution()
{
    
}


Dist_logUniform::~Dist_logUniform()
{
    
}



Dist_logUniform* Dist_logUniform::clone( void ) const
{
    return new Dist_logUniform(*this);
}


RevBayesCore::LogUniformDistribution* Dist_logUniform::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* mi   = static_cast<const Real &>( min->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* ma   = static_cast<const Real &>( max->getRevObject() ).getDagNode();
    RevBayesCore::LogUniformDistribution* d  = new RevBayesCore::LogUniformDistribution( mi, ma );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_logUniform::getClassType(void)
{
    
    static std::string rev_type = "Dist_logUniform";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_logUniform::getClassTypeSpec(void)
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
std::string Dist_logUniform::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "loguniform";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_logUniform::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Nicolas Lartillot" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_logUniform::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A strictly positive real number x has a log-uniform distribution over interval (min,max) if its logarithm y = ln(x) has uniform distribution over interval (ln(min),ln(max)).";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_logUniform::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The log-uniform distribution is defined over strictly positive real numbers.";
    details += " Saying that x is log-uniform is equivalent to saying that y = ln(x) is uniform.";
    details += " The log-uniform distribution therefore expresses lack of information about the order of magnitude of a scale parameter: ";
    details += " if x has a log-uniform distribution, then it has equal chance to be contained by any of the intervals of the form (10^k, 10^(k+1)) within the allowed range.\n\n";
    
    details += "The density is p(x) = 1/x, which can be seen by defining x = exp(y) where y has uniform distribution and apply the change-of-variable formula.\n\n";
    
    details += "The log-uniform distribution is improper when defined over the entire positive real line.";
    details += " To always make it proper, in RevBayes, a min and a max should always be specified.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_logUniform::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# a log-uniform prior over the rate of change of a Brownian trait (or a Brownian relaxed clock)\n";
    example += "trueTree = readTrees(\"data/primates.tree\")[1]\n";
    example += "sigma ~ dnLogUniform(min=0.001, max=1000)\n";
    example += "X ~ dnBrownian(trueTree,sigma)\n";
    example += "# ...\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_logUniform::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_logUniform::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnUniform" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_logUniform::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Log-Uniform Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_logUniform::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        
        dist_member_rules.push_back( new ArgumentRule( "min", RealPos::getClassTypeSpec(), "The lower bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "max", RealPos::getClassTypeSpec(), "The upper bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_logUniform::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_logUniform::printValue(std::ostream& o) const
{
    
    o << " logUniform(min=";
    if ( min != NULL ) 
    {
        o << min->getName();
    } 
    else 
    {
        o << "?";
    }
    
    o << ", max=";
    if ( max != NULL ) 
    {
        o << max->getName();
    } 
    else 
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_logUniform::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "max" )
    {
        max = var;
    }
    else if ( name == "min" )
    {
        min = var;
    }
    else 
    {
        Distribution::setConstParameter(name, var);
    }
}
