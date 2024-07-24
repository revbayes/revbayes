#include <iosfwd>
#include <string>
#include <vector>

#include "AbstractRateMatrix.h"
#include "GeneGainLossRateMatrixFunction.h"
#include "Func_GeneGainLoss.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
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
Func_GeneGainLoss::Func_GeneGainLoss( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_GeneGainLoss* Func_GeneGainLoss::clone( void ) const
{
    
    return new Func_GeneGainLoss( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_GeneGainLoss::createFunction( void ) const
{
    
    
    RevBayesCore::TypedDagNode< long >* n       = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* b     = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* d     = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    std::string m                               = static_cast<const RlString &>( this->args[3].getVariable()->getRevObject() ).getValue();

    RevBayesCore::AbstractRateMatrix::METHOD meth;
    if ( m == "scalingAndSquaring" )
    {
        meth = RevBayesCore::AbstractRateMatrix::SCALING_AND_SQUARING;
    }
    if ( m == "scalingAndSquaringPade" )
    {
        meth = RevBayesCore::AbstractRateMatrix::SCALING_AND_SQUARING_PADE;
    }
    if ( m == "scalingAndSquaringTaylor" )
    {
        meth = RevBayesCore::AbstractRateMatrix::SCALING_AND_SQUARING_TAYLOR;
    }
    if ( m == "uniformization" )
    {
        meth = RevBayesCore::AbstractRateMatrix::UNIFORMIZATION;
    }
    if ( m == "eigen" )
    {
        meth = RevBayesCore::AbstractRateMatrix::EIGEN;
    }
        
    RevBayesCore::GeneGainLossRateMatrixFunction* f = new RevBayesCore::GeneGainLossRateMatrixFunction( n->getValue(), b, d, meth );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_GeneGainLoss::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
      
        argumentRules.push_back( new ArgumentRule( "max"            , Natural::getClassTypeSpec(), "Maximum number of GeneGainLoss.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "birth"          , RealPos::getClassTypeSpec(), "Rate of gain of a single gene.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "death"          , RealPos::getClassTypeSpec(), "Rate of loss of a single gene.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        
        std::vector<std::string> optionsMethod;
        optionsMethod.push_back( "scalingAndSquaring" );
        optionsMethod.push_back( "scalingAndSquaringPade" );
        optionsMethod.push_back( "scalingAndSquaringTaylor" );
        optionsMethod.push_back( "uniformization" );
        optionsMethod.push_back( "eigen" );
        argumentRules.push_back( new OptionRule( "method", new RlString("eigen"), optionsMethod, "The method used to compute the matrix exponential." ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_GeneGainLoss::getClassType(void)
{
    
    static std::string rev_type = "Func_GeneGainLoss";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_GeneGainLoss::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_GeneGainLoss::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnGeneGainLoss";
    
    return f_name;
}


const TypeSpec& Func_GeneGainLoss::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
