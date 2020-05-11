#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_combineCharacter.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "TypeSpec.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"

using namespace RevLanguage;

/** Default constructor */
Func_combineCharacter::Func_combineCharacter( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_combineCharacter* Func_combineCharacter::clone( void ) const
{
    
    return new Func_combineCharacter( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_combineCharacter::execute( void )
{
    const AbstractHomologousDiscreteCharacterData& a = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[0].getVariable()->getRevObject() );
    const AbstractHomologousDiscreteCharacterData& b = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[1].getVariable()->getRevObject() );
        
//    size_t n = 2;
//    RevBayesCore::AbstractHomologousDiscreteCharacterData *trans_data = a.getValue().expandCharacters( n );
    RevBayesCore::AbstractHomologousDiscreteCharacterData *trans_data = a.getValue().combineCharacters( b.getValue() );
    
    return new RevVariable( new AbstractHomologousDiscreteCharacterData(trans_data) );
}


/** Get argument rules */
const ArgumentRules& Func_combineCharacter::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argumentRules.push_back( new ArgumentRule( "a", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "First character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "b", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "Second character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_combineCharacter::getClassType(void)
{
    
    static std::string rev_type = "Func_combineCharacter";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_combineCharacter::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Procedure::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_combineCharacter::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "combineCharacter";
    
    return f_name;
}

/** Get type spec */
const TypeSpec& Func_combineCharacter::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_combineCharacter::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
    
    return return_typeSpec;
}

