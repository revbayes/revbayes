#include <iosfwd>
#include <cstdint>
#include <vector>

#include "ArgumentRule.h"
#include "Func_discretizeBeta.h"
#include "Integer.h"
#include "ModelVector.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DiscretizeBetaFunction.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"


using namespace RevLanguage;

Func_discretizeBeta::Func_discretizeBeta() : TypedFunction< ModelVector< Probability > >() {
    
}

/* Clone object */
Func_discretizeBeta* Func_discretizeBeta::clone( void ) const
{
    
    return new Func_discretizeBeta( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_discretizeBeta::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* alpha = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* beta = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    bool median = static_cast<const RlBoolean &>( args[3].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::DiscretizeBetaFunction *func = new RevBayesCore::DiscretizeBetaFunction( alpha, beta, numCats, median );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeBeta::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "alpha"   , RealPos::getClassTypeSpec()  , "The alpha parameter of the beta distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"    , RealPos::getClassTypeSpec()  , "The beta parameter of the beta distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCats", Natural::getClassTypeSpec()   , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "median"  , RlBoolean::getClassTypeSpec(), "Should we use the median or mean?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeBeta::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeBeta";
    
    return rbClassName;
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeBeta::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rbClass;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeBeta::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeBeta";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeBeta::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
