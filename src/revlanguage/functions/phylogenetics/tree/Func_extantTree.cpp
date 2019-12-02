#include "Func_extantTree.h"
#include "ExtantTreeFunction.h"
#include "RlTimeTree.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_extantTree::Func_extantTree( void ) : TypedFunction<TimeTree>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_extantTree* Func_extantTree::clone( void ) const
{
    
    return new Func_extantTree( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_extantTree::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::ExtantTreeFunction* f = new RevBayesCore::ExtantTreeFunction( tau );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_extantTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_extantTree::getClassType(void)
{
    
    static std::string rev_type = "Func_extantTree";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_extantTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/**
 * Get the primary Rev name for this function.
 */
std::string Func_extantTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnExtantTree";
    
    return f_name;
}


const TypeSpec& Func_extantTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
