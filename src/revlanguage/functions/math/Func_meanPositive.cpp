#include "MeanFunction.h"
#include "Func_meanPositive.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_meanPositive::Func_meanPositive( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_meanPositive* Func_meanPositive::clone( void ) const
{
    
    return new Func_meanPositive( *this );
}


RevBayesCore::TypedFunction<double>* Func_meanPositive::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::MeanFunction* f = new RevBayesCore::MeanFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_meanPositive::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<RealPos>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_meanPositive::getClassType(void)
{
    
    static std::string rev_type = "Func_meanPositive";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_meanPositive::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_meanPositive::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "mean";
    
    return f_name;
}


const TypeSpec& Func_meanPositive::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
