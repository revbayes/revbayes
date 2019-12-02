#include "ArgumentRule.h"
#include "DiscretizeDistributionFunction.h"
#include "Ellipsis.h"
#include "Func_discretizeDistribution.h"
#include "Integer.h"
#include "ModelVector.h"
#include "RbUtil.h"
#include "Real.h"
#include "RlBoolean.h"
#include "RlContinuousDistribution.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Func_discretizeDistribution::Func_discretizeDistribution() : TypedFunction< ModelVector< Real > >()
{
    
}

/* Clone object */
Func_discretizeDistribution* Func_discretizeDistribution::clone( void ) const
{
    
    return new Func_discretizeDistribution( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_discretizeDistribution::createFunction( void ) const
{
    
    const ContinuousDistribution& rlDistribution = static_cast<const ContinuousDistribution &>( this->args[0].getVariable()->getRevObject() );
    RevBayesCore::ContinuousDistribution* g0 = static_cast<RevBayesCore::ContinuousDistribution* >( rlDistribution.createDistribution() );
    RevBayesCore::TypedDagNode<long>* numCats = static_cast<const Integer &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::DiscretizeDistributionFunction *func = new RevBayesCore::DiscretizeDistributionFunction( g0, numCats );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeDistribution::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "G0"      , ContinuousDistribution::getClassTypeSpec(), "The distribution to discretize.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "num_cats", Integer::getClassTypeSpec()               , "The number of categories into which this distribution is categorize.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeDistribution::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeDistribution";
    
    return rbClassName;
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeDistribution::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rbClass;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeDistribution::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeDistribution";
    
    return f_name;
}



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_discretizeDistribution::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_discretizeDistribution::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description;
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_discretizeDistribution::getHelpDetails(void) const
{
    // create a variable for the details of the function
    std::string details;
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_discretizeDistribution::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_discretizeDistribution::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_discretizeDistribution::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_discretizeDistribution::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "";
    
    return title;
}


/** Get type spec */
const TypeSpec& Func_discretizeDistribution::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
