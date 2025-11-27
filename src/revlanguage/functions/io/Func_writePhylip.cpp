#include <stddef.h>
#include <iosfwd>
#include <vector>
#include <ostream>

#include "ArgumentRule.h"
#include "Func_writePhylip.h"
#include "RevNullObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "RbFileManager.h"
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
Func_writePhylip* Func_writePhylip::clone( void ) const
{
    
    return new Func_writePhylip( *this );
}


/**
 * Execute the function.
 * Here we will extract the character data object from the arguments and get the file name
 * into which we shall write the character data. Then we simply create a PhylipWriter
 * instance and delegate the work
 *
 * \return NULL because the output is going into a file
 */
RevPtr<RevVariable> Func_writePhylip::execute( void )
{
    
    // get the information from the arguments for reading the file
    const std::string& file_name = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const RevBayesCore::AbstractHomologousDiscreteCharacterData &data = static_cast< const AbstractHomologousDiscreteCharacterData & >( args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::createDirectoryForFile( file_name );

    // the filestream object
    std::ofstream out_stream( file_name );
    out_stream << data.getNumberOfTaxa() << "    " << data.getNumberOfCharacters() << std::endl;
    
    const std::vector<RevBayesCore::Taxon> &taxa = data.getTaxa();
    for (std::vector<RevBayesCore::Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
    {

        if ( !data.isTaxonExcluded( it->getName() ) )
        {

            const RevBayesCore::AbstractDiscreteTaxonData &taxon = data.getTaxonData( it->getName() );

            out_stream << it->getName() << "    ";

            size_t num_chars = taxon.getNumberOfCharacters();
            for (size_t i = 0; i < num_chars; ++i)
            {
                if ( !data.isCharacterExcluded( i ) )
                {
                    const RevBayesCore::CharacterState &c = taxon.getCharacter( i );
                    out_stream << c.getStringValue();
                }
            }
            out_stream << std::endl;
        }
    }
    
    // close the stream
    out_stream.close();
    
    return NULL;
}


/**
 * Get the argument rules for this function.
 *
 * The argument rules of the writePhylip function are:
 * (1) the filename which must be a string.
 * (2) the data object that must be some character matrix.
 *
 * \return The argument rules.
 */
const ArgumentRules& Func_writePhylip::getArgumentRules( void ) const
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
const std::string& Func_writePhylip::getClassType(void)
{
    
    static std::string rev_type = "Func_writePhylip";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Func_writePhylip::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_writePhylip::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "writePhylip";
    
    return f_name;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Func_writePhylip::getTypeSpec( void ) const
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
const TypeSpec& Func_writePhylip::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
    return return_typeSpec;
}




