#include <iosfwd>
#include <cstdint>
#include <string>
#include <vector>

#include "SumIntegerFunction.h"
#include "Func_sumNatural.h"
#include "ModelVector.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_sumNatural::Func_sumNatural( void ) : TypedFunction<Natural>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_sumNatural* Func_sumNatural::clone( void ) const
{
    
    return new Func_sumNatural( *this );
}


RevBayesCore::TypedFunction<std::int64_t>* Func_sumNatural::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<std::int64_t> >* arg = static_cast<const ModelVector<Natural> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::SumIntegerFunction* f = new RevBayesCore::SumIntegerFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_sumNatural::getArgumentRules( void ) const
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


const std::string& Func_sumNatural::getClassType(void)
{
    
    static std::string rev_type = "Func_sumNatural";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_sumNatural::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<Natural>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_sumNatural::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "sum";
    
    return f_name;
}


const TypeSpec& Func_sumNatural::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
