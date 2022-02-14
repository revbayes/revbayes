#include "Func_F2x4.h"
#include "CppDoubletFuncs.h"

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
Func_F2x4* Func_F2x4::clone( void ) const
{
    return new Func_F2x4( *this );
}

typedef RevBayesCore::Simplex CSimplex;

RevBayesCore::TypedFunction< CSimplex >* Func_F2x4::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< CSimplex >* bf1 = dynamic_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< CSimplex >* bf2 = dynamic_cast<const Simplex &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr( RevBayesCore::F2x4, bf1, bf2 );
}


/* Get argument rules */
const ArgumentRules& Func_F2x4::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "baseFrequencies1", Simplex::getClassTypeSpec(), "The stationary frequencies of the 1st nucleotide.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "baseFrequencies2", Simplex::getClassTypeSpec(), "The stationary frequencies of the 2nd nucleotide.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_F2x4::getClassType(void)
{

    static std::string rev_type = "Simplex";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_F2x4::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_F2x4::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnF2x4";

    return f_name;
}


const TypeSpec& Func_F2x4::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
