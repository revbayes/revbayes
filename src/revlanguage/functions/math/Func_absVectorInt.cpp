#include <iosfwd>
#include <string>
#include <vector>

#include "AbsoluteValueVectorFunction.h"
#include "Func_absVectorInt.h"
#include "Real.h"
#include "Natural.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "GenericFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_absVectorInt::Func_absVectorInt( void ) : TypedFunction< ModelVector<Natural> >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_absVectorInt* Func_absVectorInt::clone( void ) const
{
    
    return new Func_absVectorInt( *this );
}


// Create a non-polymorphic function.
RevBayesCore::RbVector<std::int64_t>* absVectorInt(const RevBayesCore::RbVector<std::int64_t>& x)
{
    auto abs_x = new RevBayesCore::RbVector<std::int64_t>(x);
    for(auto& e: *abs_x)
        e = std::abs(e);
    return abs_x;
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<std::int64_t> >* Func_absVectorInt::createFunction( void ) const
{
    auto* arg = static_cast<const ModelVector<Integer> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::RbVector<std::int64_t> >( absVectorInt, arg );
}


/* Get argument rules */
const ArgumentRules& Func_absVectorInt::getArgumentRules( void ) const
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


const std::string& Func_absVectorInt::getClassType(void)
{
    
    static std::string rev_type = "Func_absVectorInt";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_absVectorInt::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<ModelVector<Natural> >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_absVectorInt::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "abs";
    
    return f_name;
}


const TypeSpec& Func_absVectorInt::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
