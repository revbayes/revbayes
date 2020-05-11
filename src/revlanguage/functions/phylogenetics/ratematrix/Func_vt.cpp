#include <iosfwd>

#include "Func_vt.h"
#include "RateMatrix_Vt.h"
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
Func_vt::Func_vt( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_vt* Func_vt::clone( void ) const
{
    
    return new Func_vt( *this );
}


RevPtr<RevVariable> Func_vt::execute()
{
    
    
    RevBayesCore::RateMatrix_Vt *rmj = new RevBayesCore::RateMatrix_Vt();
    RateMatrix* value = new RateMatrix( rmj );
    
    return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_vt::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


const std::string& Func_vt::getClassType(void)
{
    
    static std::string rev_type = "Func_vt";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_vt::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_vt::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnVT";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_vt::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RateGenerator::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_vt::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
