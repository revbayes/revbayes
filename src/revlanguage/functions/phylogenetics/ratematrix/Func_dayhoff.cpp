#include <iosfwd>

#include "Func_dayhoff.h"
#include "RateMatrix_Dayhoff.h"
#include "RlRateMatrix.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlRateGenerator.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_dayhoff::Func_dayhoff( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_dayhoff* Func_dayhoff::clone( void ) const
{
    
    return new Func_dayhoff( *this );
}


RevPtr<RevVariable> Func_dayhoff::execute()
{
    
    
    RevBayesCore::RateMatrix_Dayhoff *rmj = new RevBayesCore::RateMatrix_Dayhoff();
    RateMatrix* value = new RateMatrix( rmj );
    
    return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_dayhoff::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


const std::string& Func_dayhoff::getClassType(void)
{
    
    static std::string rev_type = "Func_dayhoff";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_dayhoff::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_dayhoff::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDayhoff";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_dayhoff::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RateGenerator::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_dayhoff::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
