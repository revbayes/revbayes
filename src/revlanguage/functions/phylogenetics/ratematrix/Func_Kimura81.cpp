#include "Func_Kimura81.h"

#include "Kimura81RateMatrixFunction.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;

/** default constructor */
Func_Kimura81::Func_Kimura81( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_Kimura81* Func_Kimura81::clone( void ) const
{
    
    return new Func_Kimura81( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_Kimura81::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< double >* k1 = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* k2 = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* bf = static_cast<const Simplex &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::Kimura81RateMatrixFunction* f = new RevBayesCore::Kimura81RateMatrixFunction( k1, k2, bf );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_Kimura81::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "kappa1" , RealPos::getClassTypeSpec(), "The ratio of transitions (A<->G, C<->T) to A<->C, G<->T transversions.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "kappa2" , RealPos::getClassTypeSpec(), "The ratio of A<->T, C<->G transversions to A<->C, G<->T transversions.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Simplex( RevBayesCore::RbVector<double>(4,0.25) ) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_Kimura81::getClassType(void)
{
    
    static std::string rev_type = "Func_Kimura81";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_Kimura81::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_Kimura81::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "fnKimura81" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_Kimura81::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnK81";
    
    return f_name;
}


const TypeSpec& Func_Kimura81::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
