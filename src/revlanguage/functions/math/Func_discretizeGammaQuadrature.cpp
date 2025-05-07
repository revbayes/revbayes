#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Func_discretizeGammaQuadrature.h"
#include "Integer.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DiscretizeGammaQuadratureFunction.h"
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

Func_discretizeGammaQuadrature::Func_discretizeGammaQuadrature() : TypedFunction< MatrixRealPos >() {
    
}

/* Clone object */
Func_discretizeGammaQuadrature* Func_discretizeGammaQuadrature::clone( void ) const
{
    
    return new Func_discretizeGammaQuadrature( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_discretizeGammaQuadrature::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* shape = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* rate = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::DiscretizeGammaQuadratureFunction *func = new RevBayesCore::DiscretizeGammaQuadratureFunction( shape, rate, numCats);
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeGammaQuadrature::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "shape"   , RealPos::getClassTypeSpec()  , "The shape parameter of the Gamma distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rate"    , RealPos::getClassTypeSpec()  , "The rate parameter (rate = 1/scale) of the Gamma distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCats", Integer::getClassTypeSpec()   , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeGammaQuadrature::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeGammaQuadrature";
    
	return rbClassName; 
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeGammaQuadrature::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rbClass; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeGammaQuadrature::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeGammaQuadrature";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeGammaQuadrature::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
