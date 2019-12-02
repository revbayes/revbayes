#include "MeanFunction.h"
#include "Func_mean.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_mean::Func_mean( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_mean* Func_mean::clone( void ) const
{
    
    return new Func_mean( *this );
}


RevBayesCore::TypedFunction<double>* Func_mean::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::MeanFunction* f = new RevBayesCore::MeanFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_mean::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Real>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_mean::getClassType(void)
{
    
    static std::string rev_type = "Func_mean";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_mean::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_mean::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "mean";
    
    return f_name;
}


const TypeSpec& Func_mean::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
