#include "Func_k80.h"
#include "K80RateMatrixFunction.h"
#include "RateMatrix_Kimura80.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_k80::Func_k80( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_k80* Func_k80::clone( void ) const
{
    
    return new Func_k80( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_k80::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< double >* ka = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::K80RateMatrixFunction* f = new RevBayesCore::K80RateMatrixFunction( ka );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_k80::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "kappa" , RealPos::getClassTypeSpec(), "The transition/transversion rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_k80::getClassType(void)
{
    
    static std::string rev_type = "Func_k80";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_k80::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_k80::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "fnKimura80" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_k80::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnK80";
    
    return f_name;
}


const TypeSpec& Func_k80::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
