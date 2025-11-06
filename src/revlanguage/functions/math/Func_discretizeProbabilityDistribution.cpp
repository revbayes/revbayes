#include <iosfwd>
#include <vector>

#include "ArgumentRule.h"
#include "DiscretizeDistributionFunction.h"
#include "Func_discretizeProbabilityDistribution.h"
#include "Integer.h"
#include "ModelVector.h"
#include "RlProbabilityContinuousDistribution.h"
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
#include "RbVector.h"
#include "RealPos.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class ContinuousDistribution; }


using namespace RevLanguage;

Func_discretizeProbabilityDistribution::Func_discretizeProbabilityDistribution() : TypedFunction< ModelVector< Probability > >()
{
    
}

/* Clone object */
Func_discretizeProbabilityDistribution* Func_discretizeProbabilityDistribution::clone( void ) const
{
    
    return new Func_discretizeProbabilityDistribution( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_discretizeProbabilityDistribution::createFunction( void ) const
{
       
    const ProbabilityContinuousDistribution& rlDistribution = static_cast<const ProbabilityContinuousDistribution&>( this->args[0].getVariable()->getRevObject() );
    RevBayesCore::ContinuousDistribution* g0 = static_cast<RevBayesCore::ContinuousDistribution* >( rlDistribution.createDistribution() );
    RevBayesCore::TypedDagNode<std::int64_t>* num_cats = static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::DiscretizeDistributionFunction *func = new RevBayesCore::DiscretizeDistributionFunction( g0, num_cats );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeProbabilityDistribution::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "G0"      , TypedDistribution<Probability>::getClassTypeSpec(), "The distribution to discretize.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "num_cats", Integer::getClassTypeSpec()               , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get class name of object */
const std::string& Func_discretizeProbabilityDistribution::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizeProbabilityDistribution";
    
    return rbClassName;
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeProbabilityDistribution::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rbClass;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeProbabilityDistribution::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeDistribution";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeProbabilityDistribution::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
