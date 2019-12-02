#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "Ellipsis.h"
#include "Func_readAtlas.h"
#include "ModelVector.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RevNullObject.h"
#include "RlAtlas.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlUserInterface.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "TimeAtlas.h"
#include "TimeAtlasDataReader.h"

#include <map>
#include <set>
#include <sstream>


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readAtlas* Func_readAtlas::clone( void ) const
{
    
    return new Func_readAtlas( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readAtlas::execute( void )
{
    
    // get the information from the arguments for reading the file
    const RlString& fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    
    RevBayesCore::TimeAtlasDataReader* tadr = new RevBayesCore::TimeAtlasDataReader( fn.getValue() );
    RevBayesCore::TimeAtlas* atlas = new RevBayesCore::TimeAtlas(tadr);

    return new RevVariable( new RlAtlas(atlas) );
}


/** Get argument rules */
const ArgumentRules& Func_readAtlas::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
    
        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "The name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readAtlas::getClassType(void)
{
    
    static std::string rev_type = "Func_readAtlas";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readAtlas::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readAtlas::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readAtlas";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readAtlas::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readAtlas::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlAtlas::getClassTypeSpec();
    return return_typeSpec;
}




