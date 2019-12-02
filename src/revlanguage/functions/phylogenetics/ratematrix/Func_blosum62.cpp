#include "DeterministicNode.h"
#include "Func_blosum62.h"
#include "Natural.h"
#include "RateMatrix_Blosum62.h"
#include "Real.h"
#include "RealPos.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_blosum62::Func_blosum62( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_blosum62* Func_blosum62::clone( void ) const
{
    
    return new Func_blosum62( *this );
}


RevPtr<RevVariable> Func_blosum62::execute()
{
    
    
    RevBayesCore::RateMatrix_Blosum62 *rmj = new RevBayesCore::RateMatrix_Blosum62();
    RateMatrix* value = new RateMatrix( rmj );
    
    return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_blosum62::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


const std::string& Func_blosum62::getClassType(void)
{
    
    static std::string rev_type = "Func_blosum62";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_blosum62::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/* Get return type */
const TypeSpec& Func_blosum62::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RateGenerator::getClassTypeSpec();
    
    return return_typeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_blosum62::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnBlosum62";
    
    return f_name;
}


const TypeSpec& Func_blosum62::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
