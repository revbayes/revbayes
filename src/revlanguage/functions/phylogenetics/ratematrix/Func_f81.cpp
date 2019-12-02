#include "F81RateMatrixFunction.h"
#include "Func_f81.h"
#include "RateMatrix_F81.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_f81::Func_f81( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_f81* Func_f81::clone( void ) const
{
    
    return new Func_f81( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_f81::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* bf = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::F81RateMatrixFunction* f = new RevBayesCore::F81RateMatrixFunction( bf );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_f81::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_f81::getClassType(void)
{
    
    static std::string rev_type = "Func_f81";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_f81::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_f81::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnF81";
    
    return f_name;
}


const TypeSpec& Func_f81::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
