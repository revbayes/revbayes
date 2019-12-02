#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DirichletDistribution.h"
#include "Dist_dirichlet.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlSimplex.h"
#include "StochasticNode.h"

using namespace RevLanguage;

Dist_dirichlet::Dist_dirichlet() : TypedDistribution<Simplex>()
{
    
}


Dist_dirichlet::~Dist_dirichlet()
{
    
}



Dist_dirichlet* Dist_dirichlet::clone( void ) const
{
    return new Dist_dirichlet(*this);
}


RevBayesCore::DirichletDistribution* Dist_dirichlet::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<RealPos> &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::DirichletDistribution* d              = new RevBayesCore::DirichletDistribution( a );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_dirichlet::getClassType(void)
{
    
    static std::string rev_type = "Dist_dirichlet";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_dirichlet::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Simplex>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_dirichlet::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "dirichlet";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_dirichlet::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_dirichlet::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Dirichlet probability distribution on a simplex.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_dirichlet::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The Dirichlet probability distribution is the generalization of the Beta distribution.";
    details += " A random variable from a Dirichlet distribution is a simplex, i.e., a vector of probabilities that sum to 1.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_dirichlet::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# lets get a draw from a Dirichlet distribution\n";
    example += "a <- [1,1,1,1]   # we could also use rep(1,4)\n";
    example += "b ~ dnDirichlet(a)\n";
    example += "b\n";
    example += "# let check if b really sums to 1\n";
    example += "sum(b)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_dirichlet::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_dirichlet::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "simplex" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_dirichlet::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Dirichlet Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_dirichlet::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "alpha", ModelVector<RealPos>::getClassTypeSpec(), "The concentration parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_dirichlet::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_dirichlet::printValue(std::ostream& o) const
{
    
    o << " dirichlet(alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
}


/** Set a member variable */
void Dist_dirichlet::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "alpha" ) 
    {
        alpha = var;
    }
    else 
    {
        Distribution::setConstParameter(name, var);
    }
}
