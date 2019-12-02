#include "Func_t92.h"

#include "T92RateMatrixFunction.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "Probability.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_t92::Func_t92( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_t92* Func_t92::clone( void ) const
{
    
    return new Func_t92( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_t92::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< double >* ka = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* gc = static_cast<const Probability &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::T92RateMatrixFunction* f = new RevBayesCore::T92RateMatrixFunction( gc, ka );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_t92::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "kappa" , RealPos::getClassTypeSpec(), "The transition-tranversion rate ratio.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "gc",     Probability::getClassTypeSpec(), "The frequency of GC.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_t92::getClassType(void)
{
    
    static std::string rev_type = "Func_t92";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_t92::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_t92::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnT92";
    
    return f_name;
}


const TypeSpec& Func_t92::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
