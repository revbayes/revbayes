#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "SoftBoundUniformNormalDistribution.h"
#include "ContinuousStochasticNode.h"
#include "Dist_softBoundUniformNormal.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_SoftBoundUniformNormal::Dist_SoftBoundUniformNormal() : ContinuousDistribution()
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
RevBayesCore::SoftBoundUniformNormalDistribution* Dist_SoftBoundUniformNormal::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* mi          = static_cast<const Real &>( min->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* ma          = static_cast<const Real &>( max->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s           = static_cast<const RealPos &>( sd->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* p           = static_cast<const Probability &>( prob->getRevObject() ).getDagNode();
    RevBayesCore::SoftBoundUniformNormalDistribution*   d    = new RevBayesCore::SoftBoundUniformNormalDistribution(mi, ma, s, p);
    
    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_SoftBoundUniformNormal* Dist_SoftBoundUniformNormal::clone( void ) const
{
    
    return new Dist_SoftBoundUniformNormal(*this);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_SoftBoundUniformNormal::getClassType(void)
{
    
    static std::string rev_type = "Dist_SoftBoundUniformNormal";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_SoftBoundUniformNormal::getClassTypeSpec(void)
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
std::string Dist_SoftBoundUniformNormal::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "SoftBoundUniformNormal";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_SoftBoundUniformNormal::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_SoftBoundUniformNormal::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A softbound uniform distribution with normally distributed tails outside the interval of the uniform distribution.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_SoftBoundUniformNormal::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The center piece of this distribution a uniform distribution defined between the given interval.";
    details += " A variable is drawn from that uniform distribution with probability p and with probability 1-p outside the interval.";
    details += " The probability density outside the interval is described by a normal distribution with standard deviation sd.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_SoftBoundUniformNormal::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "p ~ dnBeta(1.0,1.0)\n";
    example += "x ~ dnBernoulli(p)\n";
    example += "x.clamp(1)\n";
    example += "moves[1] = mvSlide(p, delta=0.1, weight=1.0)\n";
    example += "monitors[1] = screenmonitor(printgen=1000, separator = \"\t\", speciation)\n";
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
std::vector<RevBayesCore::RbHelpReference> Dist_SoftBoundUniformNormal::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_SoftBoundUniformNormal::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnUniform" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_SoftBoundUniformNormal::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Softbound Uniform Distribution with Normal distributed tails.";
    
    return title;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the branch rate jump process are:
 * (1) the min of the uniform distribution.
 * (2) the max of the uniform distribution.
 * (3) the standard deviation of the normal.
 * (4) the probability of belonging to distribution one.
 *
 * \return The member rules.
 */
const MemberRules& Dist_SoftBoundUniformNormal::getParameterRules(void) const
{
    
    static MemberRules distNormMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        distNormMemberRules.push_back( new ArgumentRule( "min", Real::getClassTypeSpec()       , "The min value of the uniform distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        distNormMemberRules.push_back( new ArgumentRule( "max", Real::getClassTypeSpec()       , "The max value of the uniform distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        distNormMemberRules.push_back( new ArgumentRule( "sd" , RealPos::getClassTypeSpec()    , "The standard deviation of the normal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        distNormMemberRules.push_back( new ArgumentRule( "p"  , Probability::getClassTypeSpec(), "The probability of being within the uniform distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distNormMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_SoftBoundUniformNormal::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_SoftBoundUniformNormal::printValue(std::ostream& o) const {
    
    o << "unifSBN(min=";
    if ( min != NULL ) {
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
    o << ", sd=";
    if ( sd != NULL )
    {
        o << sd->getName();
    }
    else
    {
        o << "?";
    }
    o << ", p=";
    if ( prob != NULL )
    {
        o << prob->getName();
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
void Dist_SoftBoundUniformNormal::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "min" )
    {
        min = var;
    }
    else if ( name == "max" )
    {
        max = var;
    }
    else if ( name == "sd" )
    {
        sd = var;
    }
    else if ( name == "p" )
    {
        prob = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
