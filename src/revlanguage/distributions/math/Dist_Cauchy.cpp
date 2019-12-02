#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_Cauchy.h"
#include "CauchyDistribution.h"
#include "Real.h"
#include "RealPos.h"

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_Cauchy::Dist_Cauchy() : ContinuousDistribution(){
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_Cauchy* Dist_Cauchy::clone( void ) const
{
    
    return new Dist_Cauchy(*this);
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
RevBayesCore::CauchyDistribution* Dist_Cauchy::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode<double>* l = static_cast<const Real &>( location->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s = static_cast<const RealPos &>( scale->getRevObject() ).getDagNode();
    RevBayesCore::CauchyDistribution*   d = new RevBayesCore::CauchyDistribution(l, s);
    
    return d;
}


/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Dist_Cauchy::getClassType(void)
{ 
    
    static std::string rev_type = "Dist_Cauchy";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_Cauchy::getClassTypeSpec(void)
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
 std::vector<std::string> Dist_Cauchy::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "Lorentz" );
    
    return a_names;
}

/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_Cauchy::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Cauchy";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_Cauchy::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Andrew Magee" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_Cauchy::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Cauchy distribution with location equal to ‘location’ and scale equal to ‘scale’.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_Cauchy::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The Cauchy distribution has density:\n\n";
    details += "f(x) = 1/(pi * scale) * 1/(1 + x-(location/scale)^2)";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_Cauchy::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";

    example += "# we simulate some obversations\n";
    example += "x <- rCauchy(n=10,location=0,scale=1)\n";
    example += "# let's see what the mean and the variance are.\n";
    example += "The mean will not converge with more samples, the Cauchy family has no moments.\n";
    example += "mean(x)\n";
    example += "var(x)\n";
    example += "sd(x)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_Cauchy::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
/*std::vector<std::string> Dist_Cauchy::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}
*/

/**
 * Get the title of this help entry
 */
std::string Dist_Cauchy::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Cauchy Distribution";
    
    return title;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Cauchy distribution are:
 * (1) the mean of the distribution.
 * (2) the standard deviation.
 *
 * \return The member rules.
 */
const MemberRules& Dist_Cauchy::getParameterRules(void) const
{
    
    static MemberRules distCauchyMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        distCauchyMemberRules.push_back( new ArgumentRule( "location", Real::getClassTypeSpec()   , "The location parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distCauchyMemberRules.push_back( new ArgumentRule( "scale"  , RealPos::getClassTypeSpec(), "The scale parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    
        rules_set = true;
    }
    
    return distCauchyMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_Cauchy::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_Cauchy::printValue(std::ostream& o) const
{
    
    o << " Cauchy(location=";
    if ( location != NULL )
    {
        o << location->getName();
    }
    else
    {
        o << "?";
    }
    o << ", scale=";
    if ( scale != NULL )
    {
        o << scale->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
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
void Dist_Cauchy::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "location" )
    {
        location = var;
    }
    else if ( name == "scale" )
    {
        scale = var;
    }
    else 
    {
        ContinuousDistribution::setConstParameter(name, var);
    }
}
