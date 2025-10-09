#include "Func_writeDelimitedCharacterData.h"

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "ArgumentRule.h"
#include "Delimiter.h"
#include "RevNullObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlContinuousCharacterData.h"
#include "RlString.h"
#include "DelimitedCharacterDataWriter.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ContinuousCharacterData.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

namespace RevBayesCore { class HomologousCharacterData; }



using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Func_writeDelimitedCharacterData* Func_writeDelimitedCharacterData::clone( void ) const
{
    
    return new Func_writeDelimitedCharacterData( *this );
}


/**
 * Execute the function.
 * Here we will extract the character data object from the arguments and get the file name
 * into which we shall write the character data. Then we simply create a FastaWriter
 * instance and delegate the work
 *
 * \return NULL because the output is going into a file
 */
RevPtr<RevVariable> Func_writeDelimitedCharacterData::execute( void )
{
    
    // get the information from the arguments for reading the file
    const RlString& fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    const RevBayesCore::HomologousCharacterData *data = NULL;
    
    if ( args[1].getVariable()->getRevObject().isType( AbstractHomologousDiscreteCharacterData::getClassTypeSpec() ) )
    {
        data = &static_cast< const AbstractHomologousDiscreteCharacterData & >( args[1].getVariable()->getRevObject() ).getValue();
    }
    else if ( args[1].getVariable()->getRevObject().isType( ContinuousCharacterData::getClassTypeSpec() ) )
    {
        data = &static_cast< const ContinuousCharacterData & >( args[1].getVariable()->getRevObject() ).getValue();
    }
    
    const std::string& del = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::DelimitedCharacterDataWriter writer;
    writer.writeData(fn.getValue(), *data, del);
    
    return NULL;
}


/**
 * Get the argument rules for this function.
 *
 * The argument rules of the writeCharacterDataDelimited function are:
 * (1) the filename which must be a string.
 * (2) the data object that must be some character matrix.
 *
 * \return The argument rules.
 */
const ArgumentRules& Func_writeDelimitedCharacterData::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        argument_rules.push_back( new ArgumentRule( "filename", RlString::getClassTypeSpec(), "The name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> data_arg_types;
        data_arg_types.push_back( AbstractHomologousDiscreteCharacterData::getClassTypeSpec() );
        data_arg_types.push_back( ContinuousCharacterData::getClassTypeSpec() );
        argument_rules.push_back( new ArgumentRule( "data"    , data_arg_types, "The character data object.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new Delimiter() );
        rules_set = true;
    }
    
    return argument_rules;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Func_writeDelimitedCharacterData::getClassType(void)
{
    
    static std::string rev_type = "Func_writeDelimitedCharacterData";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Func_writeDelimitedCharacterData::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_writeDelimitedCharacterData::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "writeDelimitedCharacterData";
    
    return f_name;
}


/**
 * Get the primary Rev name for this function.
 */
std::vector<std::string> Func_writeDelimitedCharacterData::getFunctionNameAliases( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::vector<std::string> f_names;
    f_names.push_back("writeCharacterDataDelimited");

    return f_names;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Func_writeDelimitedCharacterData::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Get the return type of the function.
 * This function does not return anything so the return type is NULL.
 *
 * \return NULL
 */
const TypeSpec& Func_writeDelimitedCharacterData::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
    return return_typeSpec;
}
