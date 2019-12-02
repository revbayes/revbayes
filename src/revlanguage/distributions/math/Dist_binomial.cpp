#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BinomialDistribution.h"
#include "ContinuousStochasticNode.h"
#include "Dist_binomial.h"
#include "Probability.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_binomial::Dist_binomial(void) : TypedDistribution<Natural>()
{
    
}


Dist_binomial::~Dist_binomial(void)
{
    
}



Dist_binomial* Dist_binomial::clone( void ) const
{
    
    return new Dist_binomial(*this);
}


RevBayesCore::BinomialDistribution* Dist_binomial::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<long>*    vn = static_cast<const Natural     &>( n->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* vp = static_cast<const Probability &>( p->getRevObject() ).getDagNode();
    RevBayesCore::BinomialDistribution* d  = new RevBayesCore::BinomialDistribution( vn, vp );
    return d;
}



/* Get Rev type of object */
const std::string& Dist_binomial::getClassType(void)
{
    
    static std::string rev_type = "Dist_binomial";
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_binomial::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
	return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_binomial::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "binom" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_binomial::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "binomial";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_binomial::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_binomial::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Binomial probability distribution of x successes in n trials.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_binomial::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The binomial probability distribution defines the number of success in n trials,";
    details += " where each trial has the same success probability p. The probability is given by";
    details += " (n choose x) p^(x) * (1-p)^(n-p)";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_binomial::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "p ~ dnBeta(1.0,1.0)\n";
    example += "x ~ dnBinomial(n=10,p)\n";
    example += "x.clamp(8)\n";
    example += "moves[1] = mvSlide(p, delta=0.1, weight=1.0)\n";
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
std::vector<RevBayesCore::RbHelpReference> Dist_binomial::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_binomial::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnBernoulli" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_binomial::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Binomial Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_binomial::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        dist_member_rules.push_back( new ArgumentRule( "p", Probability::getClassTypeSpec(), "Probability of success.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "n", Natural::getClassTypeSpec()    , "Number of trials.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_binomial::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_binomial::printValue(std::ostream& o) const
{
    
    o << "Binomial(p=";
    if ( p != NULL )
    {
        o << p->getName();
    }
    else
    {
        o << "?";
    }
    o << ", n=";
    if ( n != NULL )
        o << n->getName();
    else
        o << "?";
    o << ")";
}


/** Set a member variable */
void Dist_binomial::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
}
