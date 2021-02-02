#include "Func_earlyBurstRates.h"

#include "ModelVector.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "RlDeterministicNode.h"
#include "EarlyBurstRatesFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** default constructor */
Func_EarlyBurstRates::Func_EarlyBurstRates( void ) : TypedFunction<ModelVector<RealPos> >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_EarlyBurstRates* Func_EarlyBurstRates::clone( void ) const {
    
    return new Func_EarlyBurstRates( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_EarlyBurstRates::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau    = static_cast<const TimeTree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* sigma              = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* lambda             = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::EarlyBurstRatesFunction* f = new RevBayesCore::EarlyBurstRatesFunction( tau, sigma, lambda );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_EarlyBurstRates::getArgumentRules( void ) const
{
    
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set ) {
        
        argument_rules.push_back( new ArgumentRule( "tree",        TimeTree::getClassTypeSpec(), "The tree which will be re-scaled.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "sigma",       RealPos::getClassTypeSpec(), "The rate of drift.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "lambda",      RealPos::getClassTypeSpec(), "The rate of change through time.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        
        rules_set = true;
    }
    
    return argument_rules;
}


const std::string& Func_EarlyBurstRates::getClassType(void)
{
    
    static std::string rev_type = "Func_EarlyBurstRates";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_EarlyBurstRates::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_EarlyBurstRates::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnEarlyBurst";
    
    return f_name;
}


const TypeSpec& Func_EarlyBurstRates::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
