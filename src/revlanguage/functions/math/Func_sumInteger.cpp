#include "SumIntegerFunction.h"
#include "Func_sumInteger.h"
#include "ModelVector.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_sumInteger::Func_sumInteger( void ) : TypedFunction<Integer>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_sumInteger* Func_sumInteger::clone( void ) const
{
    
    return new Func_sumInteger( *this );
}


RevBayesCore::TypedFunction<long>* Func_sumInteger::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<long> >* arg = static_cast<const ModelVector<Integer> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::SumIntegerFunction* f = new RevBayesCore::SumIntegerFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_sumInteger::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Integer>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_sumInteger::getClassType(void)
{
    
    static std::string rev_type = "Func_sumInteger";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_sumInteger::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<Integer>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_sumInteger::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "sum";
    
    return f_name;
}


const TypeSpec& Func_sumInteger::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
