#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_multinomial.h"
#include "MultinomialDistribution.h"
#include "Probability.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlSimplex.h"

using namespace RevLanguage;

Dist_multinomial::Dist_multinomial(void) : TypedDistribution< ModelVector<Natural> >()
{
    
}


Dist_multinomial::~Dist_multinomial(void)
{
    
}



Dist_multinomial* Dist_multinomial::clone( void ) const
{
    
    return new Dist_multinomial(*this);
}


RevBayesCore::MultinomialDistribution* Dist_multinomial::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<long>*                              vn = static_cast<const Natural     &>( n->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >*          vp = static_cast<const Simplex &>( p->getRevObject() ).getDagNode();
    RevBayesCore::MultinomialDistribution*                        d  = new RevBayesCore::MultinomialDistribution( vp, vn );
    return d;
}



/* Get Rev type of object */
const std::string& Dist_multinomial::getClassType(void)
{
    
    static std::string rev_type = "Dist_multinomial";
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_multinomial::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Natural> >::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_multinomial::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "multinomial";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_multinomial::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_multinomial::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A multinomial distribution defines a probability distribution on a vector of natural numbers. It is understood as randomly picking n times from the k categories with replacement where each catefory has its own probability p[i].";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_multinomial::getHelpDetails(void) const
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
std::string Dist_multinomial::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "p <- simplex(1,1,1,1)\n";
    example += "x ~ dnMultinomial(10, p)\n";
    example += "y ~ dnDirichlet(x)\n";
    example += "y.clamp( simplex(1,2,3,4) )\n";
    example += "moves[1] = mvSlide(x, delta=0.1, weight=1.0)\n";
    example += "monitors[1] = screenmonitor(printgen=1000, separator = \"\t\", x)\n";
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
std::vector<RevBayesCore::RbHelpReference> Dist_multinomial::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_multinomial::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnDirichlet" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_multinomial::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Multinomial Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_multinomial::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "p", Simplex::getClassTypeSpec(), "The simplex of probabilities for the categories.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "n", Natural::getClassTypeSpec(), "The number of draws.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_multinomial::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_multinomial::printValue(std::ostream& o) const
{
    
    o << "multinomial(p=";
    if ( p != NULL )
        o << p->getName();
    else
        o << "?";
    o << ", n=";
    if ( n != NULL )
        o << n->getName();
    else
        o << "?";
    o << ")";
}


/** Set a member variable */
void Dist_multinomial::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "p" )
    {
        p = var;
    }
    else if ( name == "n" )
    {
        n = var;
    }
    else
    {
        TypedDistribution< ModelVector<Natural> >::setConstParameter(name, var);
    }
}
