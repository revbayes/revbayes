#include "DeterministicNode.h"
#include "Func_rtRev.h"
#include "Natural.h"
#include "RateMatrix_RtRev.h"
#include "Real.h"
#include "RealPos.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_rtRev::Func_rtRev( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_rtRev* Func_rtRev::clone( void ) const
{
    
    return new Func_rtRev( *this );
}


RevPtr<RevVariable> Func_rtRev::execute()
{
    
    
    RevBayesCore::RateMatrix_RtRev *rmj = new RevBayesCore::RateMatrix_RtRev();
    RateMatrix* value = new RateMatrix( rmj );
    
    return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_rtRev::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


const std::string& Func_rtRev::getClassType(void)
{
    
    static std::string rev_type = "Func_rtRev";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_rtRev::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_rtRev::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnRtRev";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_rtRev::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RateGenerator::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_rtRev::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
