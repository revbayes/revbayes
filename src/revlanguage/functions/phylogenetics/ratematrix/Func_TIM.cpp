#include "TimRateMatrixFunction.h"
#include "Func_TIM.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_TIM::Func_TIM( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_TIM* Func_TIM::clone( void ) const
{
    
    return new Func_TIM( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_TIM::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* er = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* bf = static_cast<const Simplex &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    if ( bf->getValue().size() != 4 )
    {
        throw RbException("The TIM rate matrix requires 4 base frequencies.");
    }
    
    if ( er->getValue().size() != 4 )
    {
        throw RbException("The TIM rate matrix requires 4 rate parameters.");
    }
    
    RevBayesCore::TimRateMatrixFunction* f = new RevBayesCore::TimRateMatrixFunction( er, bf );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_TIM::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "exchangeRates"  , Simplex::getClassTypeSpec(), "The exchangeability rates between states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_TIM::getClassType(void)
{
    
    static std::string rev_type = "Func_TIM";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_TIM::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_TIM::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnTIM";
    
    return f_name;
}


const TypeSpec& Func_TIM::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
