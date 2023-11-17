#include "Func_printSeed.h"

#include <cstddef>
#include <iostream>

#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "RevNullObject.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "RbHelpReference.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

Func_printSeed::Func_printSeed() : Procedure()
{
    
}

/* Clone object */
Func_printSeed* Func_printSeed::clone( void ) const
{
    
    return new Func_printSeed( *this );
}


/** Execute function: We rely on getValue and overloaded push_back to provide functionality */
RevPtr<RevVariable> Func_printSeed::execute( void )
{
    
    
    RevBayesCore::RandomNumberGenerator *rng = RevBayesCore::GLOBAL_RNG;
	unsigned int s = rng->getSeed();

	std::cout << "Current RNG Seed = " << s << "" << std::endl;
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_printSeed::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_printSeed::getClassType(void)
{
    
    static std::string rev_type = "Func_printSeed";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_printSeed::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_printSeed::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "printSeed";
    
    return f_name;
}



/** Get type spec */
const TypeSpec& Func_printSeed::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_printSeed::getReturnType( void ) const
{
    
    return RevNullObject::getClassTypeSpec();
}
