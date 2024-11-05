#include <iosfwd>
#include <vector>

#include "ArgumentRule.h"
#include "DiscretizeDistributionFunction.h"
#include "Func_discretizeDistribution.h"
#include "Integer.h"
#include "ModelVector.h"
#include "Real.h"
#include "RlContinuousDistribution.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class ContinuousDistribution; }


using namespace RevLanguage;

Func_discretizeDistribution::Func_discretizeDistribution() : TypedFunction< ModelVector< Real > >()
{
    
}

/* Clone object */
Func_discretizeDistribution* Func_discretizeDistribution::clone( void ) const
{
    
    return new Func_discretizeDistribution( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_discretizeDistribution::createFunction( void ) const
{
    
    const ContinuousDistribution& rlDistribution = static_cast<const ContinuousDistribution &>( this->args[0].getVariable()->getRevObject() );
    RevBayesCore::ContinuousDistribution* g0 = static_cast<RevBayesCore::ContinuousDistribution* >( rlDistribution.createDistribution() );
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::DiscretizeDistributionFunction *func = new RevBayesCore::DiscretizeDistributionFunction( g0, numCats );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeDistribution::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "G0"      , ContinuousDistribution::getClassTypeSpec(), "The distribution to discretize.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "num_cats", Integer::getClassTypeSpec()               , "The number of categories into which this distribution is categorize.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeDistribution::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeDistribution";
    
    return rbClassName;
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeDistribution::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rbClass;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeDistribution::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeDistribution";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeDistribution::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
