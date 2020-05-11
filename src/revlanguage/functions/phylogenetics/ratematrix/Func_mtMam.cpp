#include <iosfwd>

#include "Func_mtMam.h"
#include "RateMatrix_MtMam.h"
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
Func_mtMam::Func_mtMam( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_mtMam* Func_mtMam::clone( void ) const
{
    
    return new Func_mtMam( *this );
}


RevPtr<RevVariable> Func_mtMam::execute()
{
    
    
    RevBayesCore::RateMatrix_MtMam *rmj = new RevBayesCore::RateMatrix_MtMam();
    RateMatrix* value = new RateMatrix( rmj );
    
    return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_mtMam::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


const std::string& Func_mtMam::getClassType(void)
{
    
    static std::string rev_type = "Func_mtMam";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_mtMam::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_mtMam::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnMtMam";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_mtMam::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RateGenerator::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_mtMam::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
