#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Clade.h"
#include "Func_CladeSpecificHierarchicalBranchRate.h"
#include "ModelVector.h"
#include "CladeSpecificHierarchicalBranchRateFunction.h"
#include "RlClade.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** Default constructor */
Func_CladeSpecificHierarchicalBranchRate::Func_CladeSpecificHierarchicalBranchRate( void ) : TypedFunction< ModelVector<RealPos> >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_CladeSpecificHierarchicalBranchRate* Func_CladeSpecificHierarchicalBranchRate::clone( void ) const
{
    
    return new Func_CladeSpecificHierarchicalBranchRate( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_CladeSpecificHierarchicalBranchRate::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( args[0].getVariable()->getRevObject() ).getDagNode();
    
//    std::vector<RevBayesCore::Clade> &c = static_cast<const ModelVector<Clade> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    
    std::vector<RevBayesCore::Clade> c;
    c = static_cast<const ModelVector<Clade> &>( this->args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* br = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* cr = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::CladeSpecificHierarchicalBranchRateFunction* f = new RevBayesCore::CladeSpecificHierarchicalBranchRateFunction( tau, c, br, cr );
    
    return f;
}


/** Get argument rules */
const ArgumentRules& Func_CladeSpecificHierarchicalBranchRate::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree",       TimeTree::getClassTypeSpec(),             "The tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "clades",     ModelVector<Clade>::getClassTypeSpec(),   "The character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "baseRates",  ModelVector<RealPos>::getClassTypeSpec(), "The branch rates without clade multipliers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "cladeRates", ModelVector<RealPos>::getClassTypeSpec(), "The rate multipliers for the clades.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_CladeSpecificHierarchicalBranchRate::getClassType(void)
{
    
    static std::string rev_type = "Func_CladeSpecificHierarchicalBranchRate";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_CladeSpecificHierarchicalBranchRate::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_CladeSpecificHierarchicalBranchRate::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCladeSpecificHierarchicalBranchRate";
    
    return f_name;
}

/** Get type spec */
const TypeSpec& Func_CladeSpecificHierarchicalBranchRate::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
