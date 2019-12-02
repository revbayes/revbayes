#include "ArgumentRule.h"
#include "Dist_LKJ.h"
#include "Natural.h"
#include "RealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "StochasticNode.h"
#include "LKJDistribution.h"

using namespace RevLanguage;

Dist_LKJ::Dist_LKJ() : TypedDistribution<MatrixRealSymmetric>()
{
    
}


Dist_LKJ::~Dist_LKJ()
{
    
}



Dist_LKJ* Dist_LKJ::clone( void ) const
{
    return new Dist_LKJ(*this);
}


RevBayesCore::LKJDistribution* Dist_LKJ::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* et = static_cast<const RealPos&>( eta->getRevObject() ).getDagNode();
    size_t d  = size_t( static_cast<const Natural &>( dim->getRevObject() ).getValue() );
    
    RevBayesCore::LKJDistribution* f = new RevBayesCore::LKJDistribution( et, d );

    return f;

}



/* Get Rev type of object */
const std::string& Dist_LKJ::getClassType(void)
{
    static std::string rev_type = "Dist_LKJ";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_LKJ::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return rev_type_spec;
}

/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_LKJ::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    
    return a_names;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_LKJ::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael R. May" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_LKJ::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "The LKJ (Lewandowski-Kurowicka-Joe) distribution on correlation matrices with concentration parameter eta.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_LKJ::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "The LKJ distribution is uniform over positive-definite correlation matrices when eta=1.";
    details += "The probability density of a correlation matrix under the LKJ distribution is:";
    details += "f(x) = det(x)^(eta - 1)";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_LKJ::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "\n";
    
    example += "# we simulate a correlation matrix.\n";
    example += "R <- rLKJ(n=1, eta=1, dim=5)\n\n";
    example += "# let's print the simulated correlation matrix\n";
    example += "R\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_LKJ::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    RevBayesCore::RbHelpReference ref1;
    ref1.setCitation("Lewandowski D, D Kurowicka, H Joe (2009). Generating random correlation matrices based on vines and extended onion method. Journal of Multivariate Analysis, 100(9):1989-2001.");
    references.push_back( ref1 );

    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_LKJ::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnLKJPartial" );

    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_LKJ::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "LKJ Distribution";
    
    return title;
}



/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_LKJ::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "LKJ";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_LKJ::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "eta", RealPos::getClassTypeSpec(), "The parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "dim", Natural::getClassTypeSpec(), "The dimensions of the correlation matrix.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_LKJ::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_LKJ::printValue(std::ostream& o) const {
    
    o << " LKJ(";
    o << ")";
}


/** Set a member variable */
void Dist_LKJ::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "eta" )
    {
        eta = var;
    }
    else if ( name == "dim" )
    {
        dim = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
