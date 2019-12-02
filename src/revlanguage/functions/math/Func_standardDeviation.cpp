#include "StandardDeviationFunction.h"
#include "Func_standardDeviation.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_standardDeviation::Func_standardDeviation( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_standardDeviation* Func_standardDeviation::clone( void ) const {
    
    return new Func_standardDeviation( *this );
}


RevBayesCore::TypedFunction<double>* Func_standardDeviation::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::StandardDeviationFunction* f = new RevBayesCore::StandardDeviationFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_standardDeviation::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Real>::getClassTypeSpec(), "The vector of samples.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_standardDeviation::getClassType(void)
{
    
    static std::string rev_type = "Func_standardDeviation";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_standardDeviation::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_standardDeviation::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "stdev";
    
    return f_name;
}


const TypeSpec& Func_standardDeviation::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
