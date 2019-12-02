#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_poisson.h"
#include "PoissonDistribution.h"
#include "Natural.h"
#include "Probability.h"
#include "StochasticNode.h"

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_poisson::Dist_poisson() : TypedDistribution<Natural>() 
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
Dist_poisson* Dist_poisson::clone( void ) const 
{
    
    return new Dist_poisson(*this);
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
RevBayesCore::PoissonDistribution* Dist_poisson::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* rate    = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    RevBayesCore::PoissonDistribution* d        = new RevBayesCore::PoissonDistribution( rate );
    
    return d;
}



/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Dist_poisson::getClassType(void) 
{ 
    
    static std::string rev_type = "Dist_poisson";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_poisson::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Natural>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_poisson::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "pois" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_poisson::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "poisson";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_poisson::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_poisson::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A Poisson distribution defines probabilities for natural numbers. It is defined as the number of exponentially distributed events in a given interval.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_poisson::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details;
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_poisson::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "l ~ dnUniform(0.0,100.0)\n";
    example += "x ~ dnPoisson(l)\n";
    example += "x.clamp(10)\n";
    example += "moves[1] = mvSlide(l, delta=0.1, weight=1.0)\n";
    example += "monitors[1] = mnScreen(printgen=1000, separator = \"\t\", l)\n";
    example += "mymodel = model(l)\n";
    example += "mymcmc = mcmc(mymodel, monitors, moves)\n";
    example += "mymcmc.burnin(generations=20000,tuningInterval=100)\n";
    example += "mymcmc.run(generations=200000)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_poisson::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_poisson::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnGeom" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_poisson::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Poisson Distribution";
    
    return title;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Poisson distribution are:
 * (1) the rate lambda which must be a positive real between 0 and 1 (= a probability).
 *
 * \return The member rules.
 */
const MemberRules& Dist_poisson::getParameterRules(void) const 
{
    
    static MemberRules distPoisMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        distPoisMemberRules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The rate (rate = 1/mean) parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distPoisMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_poisson::getTypeSpec( void ) const 
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
void Dist_poisson::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "lambda" )
    {
        lambda = var;
    }
    else 
    {
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
    
}
