#include "ArgumentRule.h"
#include "Func_TaxonReader.h"
#include "ModelVector.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlUserInterface.h"
#include "TaxonReader.h"

#include <sstream>


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_TaxonReader* Func_TaxonReader::clone( void ) const
{
    
    return new Func_TaxonReader( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_TaxonReader::execute( void )
{
    
    // get the information from the arguments for reading the file
    const RlString& fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    char del = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue()[0];
    
    RevBayesCore::TaxonReader tr = RevBayesCore::TaxonReader(fn.getValue(), del);
    const std::vector<RevBayesCore::Taxon>& taxa = tr.getTaxa();
    
    return new RevVariable( new ModelVector<Taxon>( taxa ) );
}


/** Get argument rules */
const ArgumentRules& Func_TaxonReader::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        argumentRules.push_back( new ArgumentRule( "filename", RlString::getClassTypeSpec(), "Relative or absolute file name.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "delimiter", RlString::getClassTypeSpec(), "Delimiter between columns.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_TaxonReader::getClassType(void)
{
    
    static std::string rev_type = "Func_TaxonReader";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_TaxonReader::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_TaxonReader::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readTaxonData";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_TaxonReader::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_TaxonReader::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = ModelVector<Taxon>::getClassTypeSpec();
    return return_typeSpec;
}




