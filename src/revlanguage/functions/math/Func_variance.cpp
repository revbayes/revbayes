#include "VarianceFunction.h"
#include "Func_variance.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_variance::Func_variance( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_variance* Func_variance::clone( void ) const
{
    
    return new Func_variance( *this );
}


RevBayesCore::TypedFunction<double>* Func_variance::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::VarianceFunction* f = new RevBayesCore::VarianceFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_variance::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_variance::getClassType(void)
{
    
    static std::string rev_type = "Func_variance";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_variance::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<RealPos>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_variance::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "var";
    
    return f_name;
}


const TypeSpec& Func_variance::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
