#include <iosfwd>
#include <cstdint>
#include <string>
#include <vector>

#include "OrderedRateMatrixFunction.h"
#include "Func_orderedRateMatrix.h"
#include "Natural.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "OptionRule.h"
#include "RateGenerator.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_orderedRateMatrix::Func_orderedRateMatrix( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_orderedRateMatrix* Func_orderedRateMatrix::clone( void ) const
{
    
    return new Func_orderedRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_orderedRateMatrix::createFunction( void ) const
{
    
    
    RevBayesCore::TypedDagNode< std::int64_t >* n                               = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* lambda                        = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* mu                            = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Boolean >* zero_state     = static_cast<const RlBoolean &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    bool rescale                                                        = static_cast<const RlBoolean &>( this->args[4].getVariable()->getRevObject() ).getDagNode()->getValue();
    const std::string method                                            = static_cast<const RlString &>( this->args[5].getVariable()->getRevObject() ).getDagNode()->getValue();

    RevBayesCore::OrderedRateMatrixFunction* f = new RevBayesCore::OrderedRateMatrixFunction( n, lambda, mu, zero_state->getValue(), rescale, method );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_orderedRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "maxState"       , Natural::getClassTypeSpec(), "The maximum state for this rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lambda"         , RealPos::getClassTypeSpec(), "The rate of gain.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        argumentRules.push_back( new ArgumentRule( "mu"             , RealPos::getClassTypeSpec(), "The rate of loss.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        argumentRules.push_back( new ArgumentRule( "allowZeroState" , RlBoolean::getClassTypeSpec(), "Can the character go into state 0.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        argumentRules.push_back( new ArgumentRule( "rescaled",         RlBoolean::getClassTypeSpec(), "Should the matrix be normalized?", ArgumentRule::BY_VALUE,              ArgumentRule::ANY, new RlBoolean(true) ) );

        std::vector<std::string> optionsMethod;
        optionsMethod.push_back( "scalingAndSquaring" );
        optionsMethod.push_back( "scalingAndSquaringPade" );
        optionsMethod.push_back( "scalingAndSquaringTaylor" );
        optionsMethod.push_back( "uniformization" );
        optionsMethod.push_back( "eigen" );
        argumentRules.push_back( new OptionRule( "matrixExponentialMethod", new RlString("eigen"), optionsMethod, "The method used to compute the matrix exponential." ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_orderedRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_orderedRateMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_orderedRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_orderedRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnOrderedRateMatrix";
    
    return f_name;
}


const TypeSpec& Func_orderedRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
