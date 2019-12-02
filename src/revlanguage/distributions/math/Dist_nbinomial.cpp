#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "NegativeBinomialDistribution.h"
#include "ContinuousStochasticNode.h"
#include "Dist_nbinomial.h"
#include "Probability.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_nbinomial::Dist_nbinomial(void) : TypedDistribution<Natural>()
{
    
}


Dist_nbinomial::~Dist_nbinomial(void)
{
    
}



Dist_nbinomial* Dist_nbinomial::clone( void ) const
{
    
    return new Dist_nbinomial(*this);
}


RevBayesCore::NegativeBinomialDistribution* Dist_nbinomial::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<long>*   vr = static_cast<const Natural     &>( r->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* vp = static_cast<const Probability &>( p->getRevObject() ).getDagNode();
    RevBayesCore::NegativeBinomialDistribution* d  = new RevBayesCore::NegativeBinomialDistribution( vr, vp );
    return d;
}



/* Get Rev type of object */
const std::string& Dist_nbinomial::getClassType(void)
{
    
    static std::string rev_type = "Dist_nbinomial";
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_nbinomial::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
	return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_nbinomial::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "nbinom" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_nbinomial::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "nbinomial";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_nbinomial::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Walker Pett" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_nbinomial::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Negative binomial probability distribution of x successes before r failures.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_nbinomial::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The negative binomial probability distribution describes the number of successes";
    details += " before r failures, where the success probability is p. The probability is given by";
    details += " (x + r - 1 choose x) p^(x) * (1-p)^r";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_nbinomial::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "p ~ dnBeta(1.0,1.0)\n";
    example += "x ~ dnNegativeBinomial(r=10,p)\n";
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
std::vector<RevBayesCore::RbHelpReference> Dist_nbinomial::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_nbinomial::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnBinomial" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_nbinomial::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Negative Binomial Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_nbinomial::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        dist_member_rules.push_back( new ArgumentRule( "r", Natural::getClassTypeSpec(), "Number of failures.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "p", Probability::getClassTypeSpec()    , "Probability of success.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_nbinomial::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_nbinomial::printValue(std::ostream& o) const
{
    
    o << "NegativeBinomial(r=";
    if ( r != NULL )
    {
        o << r->getName();
    }
    else
    {
        o << "?";
    }
    o << ", p=";
    if ( p != NULL )
        o << p->getName();
    else
        o << "?";
    o << ")";
}


/** Set a member variable */
void Dist_nbinomial::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
        
    if ( name == "p" ) 
    {
        p = var;
    }
    else if ( name == "r" )
    {
        r = var;
    }
    else
    {
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
}
