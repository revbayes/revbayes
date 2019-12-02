#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BernoulliDistribution.h"
#include "ContinuousStochasticNode.h"
#include "Dist_bernoulli.h"
#include "Probability.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_bernoulli::Dist_bernoulli() : TypedDistribution<Natural>()
{
    
    setGuiDistributionName("Bernoulli");
    setGuiDistributionToolTip("Bernoulli distribution for 0 or 1 random variables");
}


Dist_bernoulli::~Dist_bernoulli()
{
    
}



Dist_bernoulli* Dist_bernoulli::clone( void ) const
{
    
    return new Dist_bernoulli(*this);
}


RevBayesCore::BernoulliDistribution* Dist_bernoulli::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* q   = static_cast<const Probability &>( p->getRevObject() ).getDagNode();
    RevBayesCore::BernoulliDistribution* d  = new RevBayesCore::BernoulliDistribution( q );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_bernoulli::getClassType(void)
{
    
    static std::string rev_type = "Dist_bernoulli";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_bernoulli::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_bernoulli::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "bernoulli";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_bernoulli::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "John Huelsenbeck" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_bernoulli::getHelpDescription(void) const
{
    return "A Bernoulli-distributed random variable takes the value 1 with probability p and the value 0 with probability 1-p.";
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_bernoulli::getHelpDetails(void) const
{
    return "";
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_bernoulli::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "p ~ dnBeta(1.0,1.0)\n";
    example += "x ~ dnBernoulli(p)\n";
    example += "x.clamp(1)\n";
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
std::vector<RevBayesCore::RbHelpReference> Dist_bernoulli::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_bernoulli::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnBinomial" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_bernoulli::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Bernoulli Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_bernoulli::getParameterRules(void) const
{
    
    static MemberRules distBernMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        distBernMemberRules.push_back( new ArgumentRule( "p", Probability::getClassTypeSpec(), "The probability of success.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distBernMemberRules;
}


const TypeSpec& Dist_bernoulli::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_bernoulli::printValue(std::ostream& o) const
{
    
    o << "Bernoulli(p=";
    if ( p != NULL )
    {
        o << p->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_bernoulli::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
        
    if ( name == "p" ) 
    {
        p = var;
    }
    else
    {
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
}
