#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_PhyloBranchRates.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "PhyloBranchRatesFunction.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** Default constructor */
Func_PhyloBranchRates::Func_PhyloBranchRates( void ) : TypedFunction< ModelVector<RealPos> >()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PhyloBranchRates* Func_PhyloBranchRates::clone( void ) const
{
    
    return new Func_PhyloBranchRates( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_PhyloBranchRates::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* data = static_cast<const ModelVector<Real> &>( args[1].getVariable()->getRevObject() ).getDagNode();
    bool  as_log = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::PhyloBranchRatesFunction* f = new RevBayesCore::PhyloBranchRatesFunction( tau, data, as_log );
    
    return f;
}


/** Get argument rules */
const ArgumentRules& Func_PhyloBranchRates::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "nodeStates", ModelVector<Real>::getClassTypeSpec(), "The node states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "asLog", RlBoolean::getClassTypeSpec(), "Are the node states in log-space?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_PhyloBranchRates::getClassType(void)
{
    
    static std::string rev_type = "Func_PhyloBranchRates";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_PhyloBranchRates::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction< ModelVector<RealPos> >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PhyloBranchRates::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPhyloBranchRates";
    
    return f_name;
}

std::vector<std::string> Func_PhyloBranchRates::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
//    a_names.push_back( "fnPIC" );
    
    return a_names;
}

/** Get type spec */
const TypeSpec& Func_PhyloBranchRates::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


