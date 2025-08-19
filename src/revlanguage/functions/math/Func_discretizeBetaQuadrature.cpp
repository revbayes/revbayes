#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Func_discretizeBetaQuadrature.h"
#include "Integer.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DiscretizeBetaQuadratureFunction.h"
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

Func_discretizeBetaQuadrature::Func_discretizeBetaQuadrature() : TypedFunction< MatrixRealPos >() {
    
}

/* Clone object */
Func_discretizeBetaQuadrature* Func_discretizeBetaQuadrature::clone( void ) const
{
    
    return new Func_discretizeBetaQuadrature( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_discretizeBetaQuadrature::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* alpha = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* beta = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::DiscretizeBetaQuadratureFunction *func = new RevBayesCore::DiscretizeBetaQuadratureFunction( alpha, beta, numCats);
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeBetaQuadrature::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "alpha"   , RealPos::getClassTypeSpec()  , "The alpha parameter of the beta distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"    , RealPos::getClassTypeSpec()  , "The beta parameter of the beta distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCats", Integer::getClassTypeSpec()   , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeBetaQuadrature::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeBetaQuadrature";
    
	return rbClassName; 
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeBetaQuadrature::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rbClass; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeBetaQuadrature::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeBetaQuadrature";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeBetaQuadrature::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
