#include "DeterministicNode.h"
#include "Func_mtRev.h"
#include "Natural.h"
#include "RateMatrix_MtRev24.h"
#include "Real.h"
#include "RealPos.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_mtRev::Func_mtRev( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_mtRev* Func_mtRev::clone( void ) const
{
    
    return new Func_mtRev( *this );
}


RevPtr<RevVariable> Func_mtRev::execute()
{
    
    
    RevBayesCore::RateMatrix_MtRev24 *rmj = new RevBayesCore::RateMatrix_MtRev24();
    RateMatrix* value = new RateMatrix( rmj );
    
    return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_mtRev::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


const std::string& Func_mtRev::getClassType(void)
{
    
    static std::string rev_type = "Func_mtRev";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_mtRev::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_mtRev::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnMtRev";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_mtRev::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RateGenerator::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_mtRev::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
