#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_geom.h"
#include "GeometricDistribution.h"
#include "Natural.h"
#include "Probability.h"
#include "StochasticNode.h"

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_geom::Dist_geom() : TypedDistribution<Natural>() 
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
Dist_geom* Dist_geom::clone( void ) const 
{
    
    return new Dist_geom(*this);
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
RevBayesCore::GeometricDistribution* Dist_geom::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* prob    = static_cast<const Probability &>( p->getRevObject() ).getDagNode();
    RevBayesCore::GeometricDistribution* d      = new RevBayesCore::GeometricDistribution( prob );
    
    return d;
}



/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Dist_geom::getClassType(void) 
{ 
    
    static std::string rev_type = "Dist_geom";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_geom::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Natural>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_geom::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "geom" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_geom::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "geometric";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_geom::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_geom::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A geometric distribution defines a random variable on natural numbers. The geometric distribution describes the number of success until the first failure, with success probability p.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_geom::getHelpDetails(void) const
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
std::string Dist_geom::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "p ~ dnBeta(1.0,1.0)\n";
    example += "x ~ dnGeom(p)\n";
    example += "x.clamp(10)\n";
    example += "moves[1] = mvSlide(p, delta=0.1, weight=1.0)\n";
    example += "monitors[1] = screenmonitor(printgen=1000, separator = \"\t\", p)\n";
    example += "mymodel = model(p)\n";
    example += "mymcmc = mcmc(mymodel, monitors, moves)\n";
    example += "mymcmc.burnin(generations=20000,tuningInterval=100)\n";
    example += "mymcmc.run(generations=200000)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_geom::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_geom::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnPoisson" );
    see_also.push_back( "mvRandomIntegerWalk" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_geom::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Geometric Distribution";
    
    return title;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the geometric distribution are:
 * (1) the rate lambda which must be a positive real between 0 and 1 (= a probability).
 *
 * \return The member rules.
 */
const MemberRules& Dist_geom::getParameterRules(void) const 
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        dist_member_rules.push_back( new ArgumentRule( "p", Probability::getClassTypeSpec(), "The probability of success.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_geom::getTypeSpec( void ) const 
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
void Dist_geom::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "p" ) 
    {
        p = var;
    }
    else 
    {
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
    
}
