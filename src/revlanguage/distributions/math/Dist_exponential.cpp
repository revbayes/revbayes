#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_exponential.h"
#include "ExponentialDistribution.h"
#include "Real.h"
#include "RealPos.h"

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_exponential::Dist_exponential() : PositiveContinuousDistribution() 
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
Dist_exponential* Dist_exponential::clone( void ) const 
{

    return new Dist_exponential(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be 
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the 
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::ExponentialDistribution* Dist_exponential::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* l       = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    RevBayesCore::ExponentialDistribution* d    = new RevBayesCore::ExponentialDistribution( l );
    
    return d;
}



/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Dist_exponential::getClassType(void) 
{ 
    
    static std::string rev_type = "Dist_exponential";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_exponential::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( PositiveContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_exponential::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "exp" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_exponential::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "exponential";
    
    return d_name;
}



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_exponential::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael Landis" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_exponential::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Exponential distribution with rate equal to ‘lambda’.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_exponential::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string s = "";
    s = "The exponential distribution has density:\n\n";
    s += "f(x) = 1/r * exp(-lambda*x)\n\n";
    s += "where lambda is the rate parameter.";
    
    return s;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_exponential::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# we set a rate parameter\n";
    example += "rate <- 10.0";
    example += "# we create an exponentially distributed random variable\n";
    example += "x ~ dnExponential(lambda=rate)\n";
    example += "# compute the probability of the variable\n";
    example += "x.probability()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_exponential::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_exponential::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    // see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_exponential::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Exponential Distribution";
    
    return title;
}

/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the exponential distribution are:
 * (1) the rate lambda which must be a positive real.
 *
 * \return The member rules.
 */
const MemberRules& Dist_exponential::getParameterRules(void) const 
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        dist_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The rate ( rate==1/mean) parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(1.0) ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_exponential::getTypeSpec( void ) const 
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** 
 * Set a member variable.
 * 
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_exponential::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "lambda" ) 
    {
        lambda = var;
    }
    else 
    {
        PositiveContinuousDistribution::setConstParameter(name, var);
    }
    
}
