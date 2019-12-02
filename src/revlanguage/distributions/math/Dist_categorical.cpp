#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "CategoricalDistribution.h"
#include "ContinuousStochasticNode.h"
#include "Dist_categorical.h"
#include "RlSimplex.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_categorical::Dist_categorical() :
    TypedDistribution<Natural>()
{
    setGuiDistributionName("Categorical");
    setGuiDistributionToolTip("Categorical distribution for random variables on the space {1, 2, ..., k}");
}


Dist_categorical::~Dist_categorical()
{
    
}



Dist_categorical* Dist_categorical::clone( void ) const
{
    return new Dist_categorical(*this);
}


RevBayesCore::CategoricalDistribution* Dist_categorical::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >*            q    = static_cast<const Simplex&>( p->getRevObject() ).getDagNode();
    RevBayesCore::CategoricalDistribution*                          d    = new RevBayesCore::CategoricalDistribution( q );
    
    return d;
}


/* Get Rev type of class (static) */
const std::string& Dist_categorical::getClassType( void )
{
    static std::string rev_type = "Dist_categorical";
    
    return rev_type;
}


/* Get Rev type spec of class (static) */
const TypeSpec& Dist_categorical::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_categorical::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "cat" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_categorical::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "categorical";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_categorical::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Fredrik Ronquist" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_categorical::getHelpDescription(void) const
{
    std::string description = "";
    description += "The categorical distribution, sometimes referred to as the generalized Bernoulli";
    description += " distribution. It describes the probability of one of K different outcomes,";
    description += " labeled from 1 to K, with each outcome probability separately specified.";

    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_categorical::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The argument to the constructor is a simplex containing the";
    details += " probabilities of the outcomes. The outcomes are labeled from";
    details += " 1 to K, where K is the number of elements in the simplex.";
    details += " Outcome i has probability specified by component i";
    details += " in the simplex.\n\n";
    
    details += "A typical scenario where a categorical variable is used";
    details += " is in the definition of a variable drawn from a mixture.";
    details += " A vector of mixture components is set up first, and then";
    details += " a stochastic variable drawn from a categorical distribution";
    details += " is used as an index in a deterministic assignment that";
    details += " points to a component in the mixture. See example below.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_categorical::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";

    example += "# Define a stochastic variable x that is drawn from\n";
    example += "# a categorical distribution with 4 categories, each\n";
    example += "# category having the same probability, then examine\n";
    example += "# the value of x.\n";
    example += "x ~ dnCat( simplex(1,1,1,1) )\n";
    example += "x\n";
    example += "\n";
    example += "# Draw 10 values from the distribution and place them\n";
    example += "# in a vector a, then examine a.\n";
    example += "for ( i in 1:10 ) {\n";
    example += "    a[i] <- x\n";
    example += "    x.redraw()\n";
    example += "}\n";
    example += "a\n";
    example += "\n";
    example += "# Use x in defining a deterministic variable y taking\n";
    example += "# on values from a mixture of RealPos values representing\n";
    example += "# rates from a discretized scaled gamma distribution\n";
    example += "# with four categories.\n";
    example += "shape ~ dnExp( 10.0 )\n";
    example += "rates := fnDiscretizeGamma( shape, shape, 4 )\n";
    example += "y := rates[x]\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_categorical::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_categorical::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnBinomial" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_categorical::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "The Categorical Distribution";
    
    return title;
}


/** Return member rules */
const MemberRules& Dist_categorical::getParameterRules(void) const
{
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "p", Simplex::getClassTypeSpec(), "The probability for each category.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get Rev type spec of instance (dynamic) */
const TypeSpec& Dist_categorical::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}


/** Print value for user */
void Dist_categorical::printValue(std::ostream& o) const
{
    
    o << "Categorical( p = ";
    if ( p != NULL )
    {
        if ( p->getName() != "" )
            o << p->getName();
        else
            o << "?";
    }
    else
    {
        o << "?";
    }
    o << " )";
}


/** Set a member variable */
void Dist_categorical::setConstParameter( const std::string& name, const RevPtr<const RevVariable>& var )
{
    if ( name == "p" )
    {
        p = var;
    }
    else
    {
        TypedDistribution<Natural>::setConstParameter( name, var );
    }
}
