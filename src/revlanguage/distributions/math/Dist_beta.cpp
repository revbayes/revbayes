#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BetaDistribution.h"
#include "Dist_beta.h"
#include "RealPos.h"
#include "Probability.h"
#include "RlContinuousStochasticNode.h"

using namespace RevLanguage;

Dist_beta::Dist_beta() : TypedDistribution<Probability>()
{
    
    setGuiDistributionName("Beta");
    setGuiDistributionToolTip("Beta distribution for random variables on the interval [0,1]");
}


Dist_beta::~Dist_beta()
{
    
}



Dist_beta* Dist_beta::clone( void ) const
{
    
    return new Dist_beta(*this);
}


RevBayesCore::BetaDistribution* Dist_beta::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* a   = static_cast<const RealPos &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* b   = static_cast<const RealPos &>( beta->getRevObject() ).getDagNode();
    RevBayesCore::BetaDistribution* d       = new RevBayesCore::BetaDistribution(a, b);
    
    return d;
}



Probability* Dist_beta::createRandomVariable(void) const
{
    
    RevBayesCore::ContinuousDistribution* d = createDistribution();
    RevBayesCore::TypedDagNode<double>* rv  = new ContinuousStochasticNode("", d, this->clone() );
    
    return new Probability(rv);
}



/* Get Rev type of object */
const std::string& Dist_beta::getClassType(void)
{
    
    static std::string rev_type = "Dist_beta";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_beta::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Probability>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_beta::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "beta";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_beta::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_beta::getHelpDescription(void) const
{
    return "The Beta probability distribution.";
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_beta::getHelpDetails(void) const
{
    return "The Beta distribution is a probability distribution defined on 0 to 1.";
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_beta::getHelpExample(void) const
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
std::vector<RevBayesCore::RbHelpReference> Dist_beta::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_beta::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnDirichlet" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_beta::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Beta Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_beta::getParameterRules(void) const
{
    
    static MemberRules distUnifMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        distUnifMemberRules.push_back( new ArgumentRule( "alpha", RealPos::getClassTypeSpec(), "The alpha shape parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distUnifMemberRules.push_back( new ArgumentRule( "beta" , RealPos::getClassTypeSpec(), "The beta shape parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distUnifMemberRules;
}


const TypeSpec& Dist_beta::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_beta::printValue(std::ostream& o) const
{
    
    o << "beta(alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ", beta=";
    if ( beta != NULL )
    {
        o << beta->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_beta::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
        
    if ( name == "alpha" ) 
    {
        alpha = var;
    }
    else if ( name == "beta" ) 
    {
        beta = var;
    }
    else
    {
        TypedDistribution<Probability>::setConstParameter(name, var);
    }
}
