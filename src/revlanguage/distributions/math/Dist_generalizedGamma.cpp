#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_generalizedGamma.h"
#include "GeneralizedGammaDistribution.h"
#include "Real.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_generalizedGamma::Dist_generalizedGamma() : PositiveContinuousDistribution()
{

}


Dist_generalizedGamma::~Dist_generalizedGamma()
{

}



Dist_generalizedGamma* Dist_generalizedGamma::clone( void ) const
{

    return new Dist_generalizedGamma(*this);
}


RevBayesCore::GeneralizedGammaDistribution* Dist_generalizedGamma::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* pa        = static_cast<const RealPos &>( a->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* pc        = static_cast<const RealPos &>( c->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* pl        = static_cast<const RealPos &>( l->getRevObject() ).getDagNode();
    RevBayesCore::GeneralizedGammaDistribution* d = new RevBayesCore::GeneralizedGammaDistribution(pa, pc, pl);

    return d;
}



/* Get Rev type of object */
const std::string& Dist_generalizedGamma::getClassType(void)
{

    static std::string rev_type = "Dist_generalizedGamma";

	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_generalizedGamma::getClassTypeSpec(void)
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
std::string Dist_generalizedGamma::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "generalizedGamma";

    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_generalizedGamma::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Andrew Magee" );

    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_generalizedGamma::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Generalized Gamma probability distribution for positive real numbers.";

    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_generalizedGamma::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The generalized gamma distribution (Stacy 1965) contains the gamma, Weibull, and other distributions as special cases.";
    details += " The parameterization here is that of Gomes, Combes, and Dussauchoy (2008).";

    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_generalizedGamma::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";

    example += "# lets simulate, this is equivalent to a Gamma(2,4)\n";
    example += "a <- rgeneralizedgamma(1000,a=4,c=1,l=2)\n";
    example += "# we expect a mean of 2\n";
    example += "mean(a)\n";
    example += "\n";
    example += "# create a random variable\n";
    example += "x ~ dnGeneralizedGamma(a=4,c=1,l=2)\n";
    example += "x\n";

    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_generalizedGamma::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;


    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_generalizedGamma::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnGamma" );
    see_also.push_back( "dnHalfNormal" );


    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_generalizedGamma::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Generalized Gamma Distribution";

    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_generalizedGamma::getParameterRules(void) const
{

    static MemberRules distGammaMemberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        distGammaMemberRules.push_back( new ArgumentRule( "a", RealPos::getClassTypeSpec(), "The shape parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distGammaMemberRules.push_back( new ArgumentRule( "c", RealPos::getClassTypeSpec(), "A scale parameter (if 1, distribution collapses to Gamma).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distGammaMemberRules.push_back( new ArgumentRule( "l", RealPos::getClassTypeSpec(), "A scale parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return distGammaMemberRules;
}


const TypeSpec& Dist_generalizedGamma::getTypeSpec( void ) const {

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** Print value for user */
void Dist_generalizedGamma::printValue(std::ostream& o) const {

    o << "generalizedGamma(a=";
    if ( a != NULL )
    {
        o << a->getName();
    }
    else
    {
        o << "?";
    }
    o << ", c=";
    if ( c != NULL )
    {
        o << c->getName();
    }
    else
    {
        o << "?";
    }
    if ( l != NULL )
    {
        o << l->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_generalizedGamma::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {

    if ( name == "a" )
    {
        a = var;
    }
    else if ( name == "c" )
    {
        c = var;
    }
    else if ( name == "l" )
    {
        l = var;
    }
    else
    {
        PositiveContinuousDistribution::setConstParameter(name, var);
    }
}
