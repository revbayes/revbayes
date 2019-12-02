#include <sstream>
#include <vector>

#include "ArgumentRule.h"
#include "RelativeNodeAgeConstraints.h"
#include "RelativeNodeAgeConstraintsReader.h"
#include "Func_readRelativeNodeAgeConstraints.h"
#include "RlRelativeNodeAgeConstraints.h"
#include "RlString.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/** Clone object */
Func_readRelativeNodeAgeConstraints* Func_readRelativeNodeAgeConstraints::clone( void ) const
{
    
    return new Func_readRelativeNodeAgeConstraints( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readRelativeNodeAgeConstraints::execute( void )
{
    
    // get the information from the arguments for reading the file
    const RlString& fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    
    RevBayesCore::RelativeNodeAgeConstraintsReader* dmr = new RevBayesCore::RelativeNodeAgeConstraintsReader( fn.getValue(), '\t' );
    RevBayesCore::RelativeNodeAgeConstraints* dm = new RevBayesCore::RelativeNodeAgeConstraints(dmr);
    
    return new RevVariable( new RlRelativeNodeAgeConstraints(dm) );
}


/** Get argument rules */
const ArgumentRules& Func_readRelativeNodeAgeConstraints::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "Relative or absolute name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readRelativeNodeAgeConstraints::getClassType(void)
{
    
    static std::string rev_type = "Func_readRelativeNodeAgeConstraints";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readRelativeNodeAgeConstraints::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readRelativeNodeAgeConstraints::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readRelativeNodeAgeConstraints";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readRelativeNodeAgeConstraints::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readRelativeNodeAgeConstraints::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlRelativeNodeAgeConstraints::getClassTypeSpec();
    return return_typeSpec;
}

