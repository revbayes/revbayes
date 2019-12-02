#include "Func_TamuraNei.h"
#include "TamuraNeiRateMatrixFunction.h"
#include "RateMatrix_Kimura80.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_TamuraNei::Func_TamuraNei( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_TamuraNei* Func_TamuraNei::clone( void ) const
{
    
    return new Func_TamuraNei( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_TamuraNei::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< double >* k1 = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* k2 = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* bf = static_cast<const Simplex &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TamuraNeiRateMatrixFunction* f = new RevBayesCore::TamuraNeiRateMatrixFunction( k1, k2, bf );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_TamuraNei::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "kappa1" , RealPos::getClassTypeSpec(), "The first transition rate (A<->G).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "kappa2" , RealPos::getClassTypeSpec(), "The second transition rate (C<->T).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Simplex( RevBayesCore::RbVector<double>(4,0.25) ) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_TamuraNei::getClassType(void)
{
    
    static std::string rev_type = "Func_TamuraNei";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_TamuraNei::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_TamuraNei::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "fnTamuraNei" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_TamuraNei::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnTrN";
    
    return f_name;
}


const TypeSpec& Func_TamuraNei::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
