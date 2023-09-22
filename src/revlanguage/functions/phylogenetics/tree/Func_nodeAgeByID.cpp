#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "Natural.h"
#include "Func_nodeAgeByID.h"
#include "RlBoolean.h"
#include "RlTimeTree.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "NodeAgeByID.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** default constructor */
Func_nodeAgeByID::Func_nodeAgeByID( void ) : TypedFunction<RealPos>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_nodeAgeByID* Func_nodeAgeByID::clone( void ) const {
    
    return new Func_nodeAgeByID( *this );
}


RevBayesCore::TypedFunction<double>* Func_nodeAgeByID::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    const size_t nid = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue() - 1; // this matches the indexing for the core function
    bool stemAge = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();
    RevBayesCore::NodeAgeByID* f = new RevBayesCore::NodeAgeByID( tau, nid, stemAge );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_nodeAgeByID::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree" , TimeTree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "nodeID", Natural::getClassTypeSpec()   , "The node index.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "stemAge", RlBoolean::getClassTypeSpec(), "Do we want the stem age or crown age?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_nodeAgeByID::getClassType(void)
{
    
    static std::string rev_type = "Func_nodeAgeByID";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_nodeAgeByID::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_nodeAgeByID::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "nodeAgeByID";
    
    return f_name;
}


const TypeSpec& Func_nodeAgeByID::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
