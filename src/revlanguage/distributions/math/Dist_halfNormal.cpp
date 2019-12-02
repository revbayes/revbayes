#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_halfNormal.h"
#include "HalfNormalDistribution.h"
#include "Real.h"
#include "RealPos.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_halfNormal::Dist_halfNormal() : ContinuousDistribution()
{
    
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
RevBayesCore::HalfNormalDistribution* Dist_halfNormal::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>*         o  = static_cast<const Real &>( offset->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>*         s  = static_cast<const RealPos &>( sd->getRevObject() ).getDagNode();
    RevBayesCore::HalfNormalDistribution*  d  = new RevBayesCore::HalfNormalDistribution(o, s);
    
    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_halfNormal* Dist_halfNormal::clone( void ) const
{
    
    return new Dist_halfNormal(*this);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_halfNormal::getClassType(void)
{
    
    static std::string rev_type = "Dist_halfNormal";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_halfNormal::getClassTypeSpec(void)
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
std::string Dist_halfNormal::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "halfNormal";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_halfNormal::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Andrew Magee" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_halfNormal::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "half-normal (gaussian) distribution with offset equal to ‘offset’ and standard deviation equal to ‘sd’.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_halfNormal::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The half-normal distribution has density:\n\n";
    details += " f(x) = 2/(sqrt(2 pi) sigma) e^-((x - offset)^2/(2 sigma^2))";
    details += " where offset is the offset of the distribution and sigma the standard deviation.\n\n";
    details += "f(x) = 2/(sqrt(2 pi) sigma) e^-((x - offset)^2/(2 sigma^2))\n\n";
    details += "where offset is the offset of the distribution and sigma the standard deviation.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_halfNormal::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# we simulate some oversations\n";
    example += "x <- rhalfNormal(n=10,offset=0,sd=10)\n";
    example += "# let's see what the minimum is (you could do the max too)\n";
    example += "# the minimum should never be less than the offset\n";
    example += "min(x)\n";
    example += "# let's also see what the mean and the variance are\n";
    example += "mean(x)\n";
    example += "var(x)\n";
    example += "sd(x)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_halfNormal::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_halfNormal::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnNormal" );
    see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_halfNormal::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "half-Normal Distribution";
    
    return title;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the half-normal distribution are:
 * (1) the offset of the distribution.
 * (2) the standard deviation.
 *
 * \return The member rules.
 */
const MemberRules& Dist_halfNormal::getParameterRules(void) const
{
    
    static MemberRules distHalfNormMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        distHalfNormMemberRules.push_back( new ArgumentRule( "offset", Real::getClassTypeSpec()   , "The offset parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distHalfNormMemberRules.push_back( new ArgumentRule( "sd"  , RealPos::getClassTypeSpec(), "The standard deviation parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distHalfNormMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_halfNormal::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_halfNormal::printValue(std::ostream& o) const
{
    
    o << " halfNormal(offset=";
    if ( offset != NULL )
    {
        o << offset->getName();
    }
    else
    {
        o << "?";
    }
    o << ", sd=";
    if ( sd != NULL )
    {
        o << sd->getName();
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
void Dist_halfNormal::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "offset" )
    {
        offset = var;
    }
    else if ( name == "sd" )
    {
        sd = var;
    }
    else
    {
        ContinuousDistribution::setConstParameter(name, var);
    }
}
