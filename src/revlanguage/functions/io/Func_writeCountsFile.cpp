#include <stddef.h>
#include <iosfwd>
#include <vector>

#include "ArgumentRule.h"
#include "Func_writeCountsFile.h"
#include "RevNullObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "CountsFileWriter.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }



using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Func_writeCountsFile* Func_writeCountsFile::clone( void ) const
{
    
    return new Func_writeCountsFile( *this );
}


/**
 * Execute the function.
 * Here we will extract the character data object from the arguments and get the file name
 * into which we shall write the character data. Then we simply create a CountsFileWriter
 * instance and delegate the work
 *
 * \return NULL because the output is going into a file
 */
RevPtr<RevVariable> Func_writeCountsFile::execute( void )
{
    
    // get the information from the arguments for reading the file
    const RlString& fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    const RevBayesCore::AbstractHomologousDiscreteCharacterData &data = static_cast< const AbstractHomologousDiscreteCharacterData & >( args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::CountsFileWriter fw;
    fw.writeData(fn.getValue(), data);
    
    return NULL;
}


/**
 * Get the argument rules for this function.
 *
 * The argument rules of the writeCountsFile function are:
 * (1) the filename which must be a string.
 * (2) the data object that must be some character matrix.
 *
 * \return The argument rules.
 */
const ArgumentRules& Func_writeCountsFile::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        argumentRules.push_back( new ArgumentRule( "filename", RlString::getClassTypeSpec(), "The name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "data"    , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The character data object.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Func_writeCountsFile::getClassType(void)
{
    
    static std::string rev_type = "Func_writeCountsFile";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Func_writeCountsFile::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_writeCountsFile::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "writeCountsFile";
    
    return f_name;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Func_writeCountsFile::getTypeSpec( void ) const
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
const TypeSpec& Func_writeCountsFile::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
    return return_typeSpec;
}




