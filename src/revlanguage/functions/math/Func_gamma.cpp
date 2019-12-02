#include "Func_gamma.h"
#include "GammaFunction.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_gamma::Func_gamma( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_gamma* Func_gamma::clone( void ) const
{
    
    return new Func_gamma( *this );
}


RevBayesCore::TypedFunction<double>* Func_gamma::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::GammaFunction* f = new RevBayesCore::GammaFunction( x );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_gamma::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x"    , Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_gamma::getClassType(void)
{
    
    static std::string rev_type = "Func_gamma";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_gamma::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_gamma::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "gamma";
    
    return f_name;
}


const TypeSpec& Func_gamma::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
