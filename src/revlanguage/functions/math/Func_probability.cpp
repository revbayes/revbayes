#include "Func_probability.h"
#include "Probability.h"
#include "Real.h"
#include "ReferenceFunction.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_probability::Func_probability( void ) : TypedFunction<Probability>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_probability* Func_probability::clone( void ) const
{
    
    return new Func_probability( *this );
}


RevBayesCore::TypedFunction<double>* Func_probability::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::ReferenceFunction<double>* f = new RevBayesCore::ReferenceFunction<double>( x );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_probability::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_probability::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_probability::getHelpDescription(void) const
{
    std::string description = "";
    description += "This Rev function converts between Rev types from Real to Probability.";
    
    return description;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_probability::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# Create a RealPos\n";
    example += "x <- 12/13\n";
    example += "type(x)\n";
    example += "\n";
    example += "# Convert the RealPos to Probability\n";
    example += "y := Probability(x)\n";
    example += "type(y)\n";
    
    return example;
}


const std::string& Func_probability::getClassType(void)
{
    
    static std::string rev_type = "Func_probability";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_probability::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_probability::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "Probability";
    
    return f_name;
}


const TypeSpec& Func_probability::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
