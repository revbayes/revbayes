#include "Func_stochasticMatrix.h"
#include "ArgumentRule.h"
#include "ModelVector.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "StochasticMatrixFunction.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "Real.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlSimplex.h"
#include "RlStochasticMatrix.h"
#include "Simplex.h"
#include "StringUtilities.h"


using namespace RevLanguage;


/** Default constructor */
Func_stochasticMatrix::Func_stochasticMatrix( void ) : TypedFunction< StochasticMatrix >()
{
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_stochasticMatrix* Func_stochasticMatrix::clone( void ) const
{
    return new Func_stochasticMatrix( *this );
}


/** Execute function: Construct stochastic matrix from the vector of simplices. */
RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_stochasticMatrix::createFunction( void ) const
{

    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::Simplex > >* mtx = static_cast< const ModelVector< Simplex >& >( args[0].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::StochasticMatrixFunction* func = new RevBayesCore::StochasticMatrixFunction( mtx );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_stochasticMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", ModelVector< Simplex >::getClassTypeSpec(), "A vector of simplices. Each simplex corresponds to a row in the stochastic matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object (static) */
const std::string& Func_stochasticMatrix::getClassType( void )
{
    static std::string rev_type = "Func_stochasticMatrix";
    
    return rev_type;
}


/** Get Rev type spec of object (static) */
const TypeSpec& Func_stochasticMatrix::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_stochasticMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "stochasticMatrix";
    
    return f_name;
}


/** Get Rev type spec of object from an instance */
const TypeSpec& Func_stochasticMatrix::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

