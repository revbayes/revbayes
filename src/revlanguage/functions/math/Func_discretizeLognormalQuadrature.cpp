#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Func_discretizeLognormalQuadrature.h"
#include "Integer.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DiscretizeLognormalQuadratureFunction.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlMatrixRealPos.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"


using namespace RevLanguage;

Func_discretizeLognormalQuadrature::Func_discretizeLognormalQuadrature() : TypedFunction< MatrixRealPos >() {
    
}

/* Clone object */
Func_discretizeLognormalQuadrature* Func_discretizeLognormalQuadrature::clone( void ) const
{
    
    return new Func_discretizeLognormalQuadrature( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_discretizeLognormalQuadrature::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* m   = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s   = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::DiscretizeLognormalQuadratureFunction *func = new RevBayesCore::DiscretizeLognormalQuadratureFunction( m, s, numCats);
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeLognormalQuadrature::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "mean"   , Real::getClassTypeSpec()  , "The mean (in log-space) of the Lognormal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "sd"    , RealPos::getClassTypeSpec()  , "The standard deviation (in log-space) of the Lognormal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCats", Integer::getClassTypeSpec()   , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeLognormalQuadrature::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeLognormalQuadrature";
    
	return rbClassName; 
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeLognormalQuadrature::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rbClass; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeLognormalQuadrature::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeLognormalQuadrature";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeLognormalQuadrature::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
