#include "Func_F1x4.h"
#include "CppCodonFuncs.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Simplex.h"
#include "TypeSpec.h"
#include "CodonState.h"

using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_F1x4* Func_F1x4::clone( void ) const
{
    return new Func_F1x4( *this );
}

typedef RevBayesCore::Simplex CSimplex;

CSimplex* F1x4(const CSimplex& nuc_pi1)
{
    return RevBayesCore::F1x4(nuc_pi1).clone();
}

RevBayesCore::TypedFunction< CSimplex >* Func_F1x4::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< CSimplex >* bf = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    if ( bf->getValue().size() != 4 )
    {
        throw RbException("The fnF1x4 function takes 4 base frequencies.");
    }

    return RevBayesCore::generic_function_ptr< CSimplex >( F1x4, bf );
}


/* Get argument rules */
const ArgumentRules& Func_F1x4::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the nucleotides.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_F1x4::getClassType(void)
{

    static std::string rev_type = "Simplex";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_F1x4::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_F1x4::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnF1x4";

    return f_name;
}


const TypeSpec& Func_F1x4::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
