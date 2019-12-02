#include "Func_freeSymmetricRateMatrix.h"

#include "FreeSymmetricRateMatrixFunction.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_freeSymmetricRateMatrix::Func_freeSymmetricRateMatrix( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_freeSymmetricRateMatrix* Func_freeSymmetricRateMatrix::clone( void ) const
{
    
    return new Func_freeSymmetricRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_freeSymmetricRateMatrix::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* tr = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    bool r = static_cast<const RlBoolean &>( this->args[1].getVariable()->getRevObject() ).getDagNode()->getValue();
    const std::string method = static_cast<const RlString &>( this->args[2].getVariable()->getRevObject() ).getDagNode()->getValue();
    
    RevBayesCore::FreeSymmetricRateMatrixFunction* f = new RevBayesCore::FreeSymmetricRateMatrixFunction( tr, r, method);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_freeSymmetricRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "transition_rates",   ModelVector<RealPos>::getClassTypeSpec(),   "The transition rates between states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rescaled",          RlBoolean::getClassTypeSpec(),              "Should the matrix be normalized?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        std::vector<std::string> optionsMethod;
        optionsMethod.push_back( "scalingAndSquaring" );
        optionsMethod.push_back( "scalingAndSquaringPade" );
        optionsMethod.push_back( "scalingAndSquaringTaylor" );
        optionsMethod.push_back( "uniformization" );
        optionsMethod.push_back( "eigen" );
        argumentRules.push_back( new OptionRule( "matrixExponentialMethod", new RlString("scalingAndSquaring"), optionsMethod, "The method used to compute the matrix exponential." ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_freeSymmetricRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_freeSymmetricRateMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_freeSymmetricRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_freeSymmetricRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnFreeSymmetricRateMatrix";
    
    return f_name;
}


const TypeSpec& Func_freeSymmetricRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
