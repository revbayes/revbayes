#include <iosfwd>
#include <string>
#include <vector>

#include "CAFERateMatrixFunction.h"
#include "Func_CAFE.h"
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
#include "RateGenerator.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_CAFE::Func_CAFE( void ) : TypedFunction<RateMatrix>( ) 
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_CAFE* Func_CAFE::clone( void ) const
{
    
    return new Func_CAFE( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_CAFE::createFunction( void ) const 
{    
    
    long n                                      = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< double >* birth = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* death = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::CAFERateMatrixFunction* f = new RevBayesCore::CAFERateMatrixFunction( n, birth, death );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_CAFE::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
      
        argumentRules.push_back( new ArgumentRule( "max"            , Natural::getClassTypeSpec(), "Maximum number of CAFE.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "birth"          , RealPos::getClassTypeSpec(), "Rate of gain of a single gene.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "death"          , RealPos::getClassTypeSpec(), "Rate of loss of a single gene.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_CAFE::getClassType(void)
{
    
    static std::string rev_type = "Func_CAFE";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_CAFE::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_CAFE::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCAFE";
    
    return f_name;
}


const TypeSpec& Func_CAFE::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
