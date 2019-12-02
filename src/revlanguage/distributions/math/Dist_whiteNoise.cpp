#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_whiteNoise.h"
#include "Real.h"
#include "RealPos.h"

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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_whiteNoise::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_whiteNoise::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "White-Noise process for positive real numbers.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_whiteNoise::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The white-noise process is a process of a positive continuous variable similar to Brownian motion and the Ornstein-Uhlenbeck process.";
    details += " However, the white-noise process has a large variance when the time is small, and has small variance if the time is large.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_whiteNoise::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# lets simulate\n";
    example += "a <- rWhiteNoise(1000,mu=1,sigma=4,time=4)\n";
    example += "# we expect a mean of 1\n";
    example += "mean(a)\n";
    example += "\n";
    example += "# create a random variable\n";
    example += "x ~ dnWhiteNoise(mu=1.0,sigma=4,time=1)\n";
    example += "x\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_whiteNoise::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_whiteNoise::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnGamma" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_whiteNoise::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "White-Noise Process";
    
    return title;
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
