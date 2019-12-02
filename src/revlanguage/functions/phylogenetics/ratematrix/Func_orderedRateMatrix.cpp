#include "OrderedRateMatrixFunction.h"
#include "Func_orderedRateMatrix.h"
#include "Natural.h"
#include "RateMatrix_Ordered.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

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
    
    
    RevBayesCore::TypedDagNode< long >* n          = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* lambda   = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* mu       = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::OrderedRateMatrixFunction* f = new RevBayesCore::OrderedRateMatrixFunction( n, lambda, mu );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_orderedRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "maxState"   , Natural::getClassTypeSpec(), "The maximum state for this rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lambda"     , RealPos::getClassTypeSpec(), "The rate of gain.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        argumentRules.push_back( new ArgumentRule( "mu"         , RealPos::getClassTypeSpec(), "The rate of loss.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
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
