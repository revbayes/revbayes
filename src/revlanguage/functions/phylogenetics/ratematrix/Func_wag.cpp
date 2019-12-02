#include "WagRateMatrixFunction.h"
#include "Func_wag.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "Simplex.h"
#include "TypedDagNode.h"

#include <iostream>

using namespace RevLanguage;

/** default constructor */
Func_wag::Func_wag( void ) : TypedFunction<RateMatrix>( )
{
    std::vector<double> bf_empirical_init(20);
    bf_empirical_init[ 0] = 0.0866279;
    bf_empirical_init[ 1] = 0.0439720;
    bf_empirical_init[ 2] = 0.0390894;
    bf_empirical_init[ 3] = 0.0570451;
    bf_empirical_init[ 4] = 0.0193078;
    bf_empirical_init[ 5] = 0.0367281;
    bf_empirical_init[ 6] = 0.0580589;
    bf_empirical_init[ 7] = 0.0832518;
    bf_empirical_init[ 8] = 0.0244313;
    bf_empirical_init[ 9] = 0.048466;
    bf_empirical_init[10] = 0.086209;
    bf_empirical_init[11] = 0.0620286;
    bf_empirical_init[12] = 0.0195027;
    bf_empirical_init[13] = 0.0384319;
    bf_empirical_init[14] = 0.0457631;
    bf_empirical_init[15] = 0.0695179;
    bf_empirical_init[16] = 0.0610127;
    bf_empirical_init[17] = 0.0143859;
    bf_empirical_init[18] = 0.0352742;
    bf_empirical_init[19] = 0.0708956;
    
    bf_empirical = new Simplex( RevBayesCore::Simplex(bf_empirical_init));
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_wag* Func_wag::clone( void ) const
{
    
    return new Func_wag( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_wag::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* bf = NULL;
    const RevObject* obj = &this->args[0].getVariable()->getRevObject();
    if (obj == NULL) {
        std::cout << "base frequencies equal default argument\n";
        
        // populate default empirical base frequencies
        ;
    }
    else {
        bf = static_cast<const Simplex &>( *obj ).getDagNode();
    }
    
    
    
    if (bf->getValue().size() != 20)
    {
        throw RbException("Amino acid base frequencies must size 20.");
    }
    
    RevBayesCore::WagRateMatrixFunction* f = new RevBayesCore::WagRateMatrixFunction( bf );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_wag::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, bf_empirical ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_wag::getClassType(void)
{
    
    static std::string rev_type = "Func_wag";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_wag::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_wag::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnWAG";
    
    return f_name;
}


const TypeSpec& Func_wag::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
