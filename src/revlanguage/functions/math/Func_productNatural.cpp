#include <iosfwd>
#include <string>
#include <vector>

#include "ProductIntegerFunction.h"
#include "Func_productNatural.h"
#include "ModelVector.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_productNatural::Func_productNatural( void ) : TypedFunction<Natural>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_productNatural* Func_productNatural::clone( void ) const
{
    
    return new Func_productNatural( *this );
}


RevBayesCore::TypedFunction<long>* Func_productNatural::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<long> >* arg = static_cast<const ModelVector<Natural> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::ProductIntegerFunction* f = new RevBayesCore::ProductIntegerFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_productNatural::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Natural>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_productNatural::getClassType(void)
{
    
    static std::string rev_type = "Func_productNatural";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_productNatural::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_productNatural::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "prod";
    
    return f_name;
}


const TypeSpec& Func_productNatural::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
