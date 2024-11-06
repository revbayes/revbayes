#include <iosfwd>
#include <cstdint>
#include <vector>

#include "ArgumentRule.h"
#include "Func_discretizeGamma.h"
#include "Integer.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DiscretizeGammaFunction.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
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

Func_discretizeGamma::Func_discretizeGamma() : TypedFunction< ModelVector< RealPos > >() {
    
}

/* Clone object */
Func_discretizeGamma* Func_discretizeGamma::clone( void ) const
{
    
    return new Func_discretizeGamma( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_discretizeGamma::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* shape = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* rate = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    bool median = static_cast<const RlBoolean &>( args[3].getVariable()->getRevObject() ).getValue();

    RevBayesCore::DiscretizeGammaFunction *func = new RevBayesCore::DiscretizeGammaFunction( shape, rate, numCats, median );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeGamma::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "shape"   , RealPos::getClassTypeSpec()  , "The shape parameter of the gamma distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rate"    , RealPos::getClassTypeSpec()  , "The rate parameter (rate = 1/scale) of the gamma distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCats", Integer::getClassTypeSpec()   , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "median"  , RlBoolean::getClassTypeSpec(), "Should we use the median or mean?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeGamma::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeGamma";
    
	return rbClassName; 
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeGamma::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rbClass; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeGamma::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeGamma";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeGamma::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
