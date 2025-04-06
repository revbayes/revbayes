#include <iosfwd>
#include <string>
#include <vector>

#include "SegregatingSitesFunction.h"
#include "Func_SegregatingSites.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "Natural.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }

using namespace RevLanguage;

/** default constructor */
Func_SegregatingSites::Func_SegregatingSites( void ) : TypedFunction<Natural>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_SegregatingSites* Func_SegregatingSites::clone( void ) const
{
    
    return new Func_SegregatingSites( *this );
}


RevBayesCore::TypedFunction< std::int64_t >* Func_SegregatingSites::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData >* d = static_cast<const AbstractHomologousDiscreteCharacterData &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Boolean >* ex = static_cast<const RlBoolean &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::SegregatingSitesFunction* f = new RevBayesCore::SegregatingSitesFunction( d, ex->getValue() );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_SegregatingSites::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "data", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The alignment for which to compute the number of segregating sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "excludeAmbiguous", RlBoolean::getClassTypeSpec(), "Should we exclude ambiguous or missing characters?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_SegregatingSites::getClassType(void)
{
    
    static std::string rev_type = "Func_SegregatingSites";
    
    return rev_type;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_SegregatingSites::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnSegregatingSites";
    
    return f_name;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_SegregatingSites::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<Natural>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


const TypeSpec& Func_SegregatingSites::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
