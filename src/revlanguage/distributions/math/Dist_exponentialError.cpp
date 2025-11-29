#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Dist_exponentialError.h"
#include "Natural.h"
#include "RealPos.h"
#include "StochasticNode.h"
#include "AverageDistanceMatrix.h"
#include "DistanceMatrix.h"
#include "ExponentialErrorDistribution.h"
#include "ArgumentRules.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAverageDistanceMatrix.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistanceMatrix.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/**
* Default constructor.
*
* The default constructor does nothing except allocating the object.
*/
Dist_exponentialError::Dist_exponentialError() : TypedDistribution<AverageDistanceMatrix>()
{
    
}


Dist_exponentialError::~Dist_exponentialError()
{
    
}


/**
* Clone the object
*
* \return a clone of the object.
*/
Dist_exponentialError* Dist_exponentialError::clone( void ) const
{
    return new Dist_exponentialError(*this);
}


RevBayesCore::ExponentialErrorDistribution* Dist_exponentialError::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::AverageDistanceMatrix>* averDM = NULL;
    RevBayesCore::TypedDagNode<RevBayesCore::DistanceMatrix>* DM = NULL;
    RevBayesCore::TypedDagNode<double>* la = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    
    if ( avgDistMat->getRevObject() != RevNullObject::getInstance() )
    {
        averDM = static_cast<const AverageDistanceMatrix &>( avgDistMat->getRevObject() ).getDagNode();
    }
    if ( distMat->getRevObject() != RevNullObject::getInstance() )
    {
        DM = static_cast<const DistanceMatrix &>( distMat->getRevObject() ).getDagNode();
    }
    
    RevBayesCore::ExponentialErrorDistribution* ee = NULL;
    
    if ( averDM != NULL )
    {
        // parameter is avgDistMat
        ee = new RevBayesCore::ExponentialErrorDistribution( averDM, la );
    }
    else
    {
        // parameter is distMat
        ee = new RevBayesCore::ExponentialErrorDistribution( DM, la );
    }
    
    return ee;
    
}


/* Get class name of object */
const std::string& Dist_exponentialError::getClassType(void) {
    
    static std::string revClassType = "Dist_exponentialError";
    
    return revClassType;
}


/* Get class type spec describing type of object */
const TypeSpec& Dist_exponentialError::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return revClassTypeSpec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_exponentialError::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "expError" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_exponentialError::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "ExponentialError";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_exponentialError::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "David Cerny" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_exponentialError::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Exponential error distribution on distance matrices: the element-wise sum of squared differences from the reference matrix (supplied in the parameter avgDistMat or distMat) is exponentially distributed with rate lambda.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_exponentialError::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "The probability density of a distance matrix under the exponential error distribution is: ";
    details += "f(x) = \\lambda exp{-\\lambda \\sum_i^n \\sum_j^n [x(i,j) - D(i,j)]^2}, ";
    details += "where x(i,j) is the path-length distance between i and j in matrix x and D(i,j) is their path-length distance in the reference matrix.";
    details += "If the reference matrix is an average distance matrix, a Boolean mask is applied to ensure that only meaningful distances are used.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_exponentialError::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "\n";
    
    example += "# we read in the reference matrix\n";
    example += "DM <- readDistanceMatrix(...)\n\n";
    example += "# we simulate a distance matrix from the distribution\n";
    example += "d ~ dnExponentialError(DM, 1)\n\n";
    example += "# calculate the probability density of the simulated distance matrix\n";
    example += "d.lnProbability()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 */
std::vector<RevBayesCore::RbHelpReference> Dist_exponentialError::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    RevBayesCore::RbHelpReference ref1;
    ref1.setCitation("Lapointe F-J, C Levasseur (2004). Everything you always wanted to know about the average consensus, and more. 87-105 in Bininda-Emonds ORP, ed. Phylogenetic Supertrees: Combining Information to Reveal the Tree of Life. Computational Biology, vol. 3. Dordrecht: Kluwer Academic Publishers.");
    references.push_back( ref1 );
    
    RevBayesCore::RbHelpReference ref2;
    ref2.setCitation("Steel M, A Rodrigo (2008). Maximum likelihood supertrees. Systematic Biology, 57(2):243-250.");
    references.push_back( ref2 );

    return references;
}


/**
 * Get the title of this help entry
 */
std::string Dist_exponentialError::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Exponential Error Distribution";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_exponentialError::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "avgDistMat", AverageDistanceMatrix::getClassTypeSpec(), "The average distance matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "distMat", DistanceMatrix::getClassTypeSpec(), "The distance matrix of the proposed supertree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The rate of the exponential distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_exponentialError::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_exponentialError::printValue(std::ostream& o) const
{
    
    o << " ExponentialError(";
    o << ")";
}


/** Set a member variable */
void Dist_exponentialError::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "avgDistMat" )
    {
        avgDistMat = var;
    }
    else if ( name == "distMat" )
    {
        distMat = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
