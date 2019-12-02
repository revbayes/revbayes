#include "Func_hky.h"
#include "HkyRateMatrixFunction.h"
#include "RateMatrix_HKY.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_hky::Func_hky( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_hky* Func_hky::clone( void ) const
{
    
    return new Func_hky( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_hky::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< double >* ka = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* bf = static_cast<const Simplex &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    if ( bf->getValue().size() != 4 )
    {
        throw RbException("The HKY rate matrix can only use 4 base fequencies.");
    }
    
    RevBayesCore::HkyRateMatrixFunction* f = new RevBayesCore::HkyRateMatrixFunction( ka, bf );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_hky::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "kappa"          , RealPos::getClassTypeSpec(), "The transition-transversion rate ratio.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_hky::getClassType(void)
{
    
    static std::string rev_type = "Func_hky";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_hky::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_hky::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnHKY";
    
    return f_name;
}


const TypeSpec& Func_hky::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
