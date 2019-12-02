#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_lnorm.h"
#include "LognormalDistribution.h"
#include "Real.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_lnorm::Dist_lnorm() : PositiveContinuousDistribution() {
    
}


Dist_lnorm* Dist_lnorm::clone( void ) const {
    return new Dist_lnorm(*this);
}


RevBayesCore::LognormalDistribution* Dist_lnorm::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode<double>* m   = static_cast<const Real &>( mean->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s   = static_cast<const RealPos &>( sd->getRevObject() ).getDagNode();
    RevBayesCore::LognormalDistribution* d  = new RevBayesCore::LognormalDistribution(m, s);
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_lnorm::getClassType(void) { 
    
    static std::string rev_type = "Dist_lnorm";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_lnorm::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( PositiveContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_lnorm::getDistributionFunctionAliases( void ) const
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
std::string Dist_lnorm::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "lognormal";
    
    return d_name;
}



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_lnorm::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael Landis" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_lnorm::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Lognormal distribution is the distribution for a log-transformed normally distributed random variable with mean 'mu' and standard deviation 'sigma'.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_lnorm::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The lognormal random variable is defined as\n\n:";
    
    details += "X = exp(mu + sigma Z)\n\n";
    
    details += "where mu is the mean parameter, sigma is the standard deviation, and Z is a standard normal random variable. ";
    details += "Note, in effect, the mean and standard deviation provide the location and scale of the exponentiated normal variate, mu + sigma Z.";
    details += "The lognormal distribution has density:\n\n";
    
    details += "f(x) = 1/(x sigma sqrt(2 pi)) e^-((ln x - mu)^2/(2 sigma^2))\n\n";
    
    details += "where mu is the mean of the distribution and sigma the standard deviation.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_lnorm::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# set an expected value for x\n";
    example += "expectation_of_x <- 1\n";
    example += "# set a mean and sd parameter\n";
    example += "sd <- 0.5\n";
    example += "mean <- ln(expectation_of_x) - 0.5 * sd^2\n";
    example += "# create a lognormal distribution with expected value of 1\n";
    example += "x ~ dnLognormal(mean=mean, sd=sd)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_lnorm::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_lnorm::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    // see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_lnorm::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Lognormal Distribution";
    
    return title;
}

/** Return member rules (no members) */
const MemberRules& Dist_lnorm::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        
        memberRules.push_back( new ArgumentRule( "mean",   Real::getClassTypeSpec()   , "The mean in log-space (observed mean is exp(m))." , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "sd"  ,   RealPos::getClassTypeSpec(), "The standard deviation in log-space."             , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_lnorm::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_lnorm::printValue(std::ostream& o) const {
    
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
    o << ")";
}


/** Set a member variable */
void Dist_lnorm::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "mean" ) 
    {
        mean = var;
    }
    else if ( name == "sd" ) 
    {
        sd = var;
    }
    else 
    {
        PositiveContinuousDistribution::setConstParameter(name, var);
    }
}
