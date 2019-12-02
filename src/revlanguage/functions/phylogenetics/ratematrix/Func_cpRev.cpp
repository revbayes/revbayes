#include "DeterministicNode.h"
#include "Func_cpRev.h"
#include "Natural.h"
#include "RateMatrix_CpRev.h"
#include "Real.h"
#include "RealPos.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_cpRev::Func_cpRev( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_cpRev* Func_cpRev::clone( void ) const
{
    
    return new Func_cpRev( *this );
}


RevPtr<RevVariable> Func_cpRev::execute()
{
    
    
    RevBayesCore::RateMatrix_CpRev *rmj = new RevBayesCore::RateMatrix_CpRev();
    RateMatrix* value = new RateMatrix( rmj );
    
    return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_cpRev::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


const std::string& Func_cpRev::getClassType(void)
{
    
    static std::string rev_type = "Func_cpRev";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_cpRev::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/* Get return type */
const TypeSpec& Func_cpRev::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RateGenerator::getClassTypeSpec();
    
    return return_typeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_cpRev::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCpRev";
    
    return f_name;
}


const TypeSpec& Func_cpRev::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
