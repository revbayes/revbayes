#include <string>
#include <iosfwd>
#include <vector>

#include "Func_MinBLTimeScaling.h"

#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "ModelVector.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "TypeSpec.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "RlTree.h"

using namespace RevLanguage;

/** Default constructor */
Func_MinBLTimeScaling::Func_MinBLTimeScaling( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objects.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_MinBLTimeScaling* Func_MinBLTimeScaling::clone( void ) const
{
    
    return new Func_MinBLTimeScaling( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_MinBLTimeScaling::execute( void )
{
    RevBayesCore::Tree* treeToScale              = static_cast<const Tree&>( this->args[0].getVariable()->getRevObject() ).getValue().clone();
    std::vector<RevBayesCore::Taxon> taxonVector = static_cast<const ModelVector<Taxon> &>( this->args[1].getVariable()->getRevObject() ).getValue();
    double minBranchLength                       = static_cast<const RealPos&>( this->args[2].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::Tree *my_tree = RevBayesCore::TreeUtilities::minBLTimeScaling( *treeToScale, taxonVector, minBranchLength );
    
    return new RevVariable( new TimeTree( my_tree ) );
}


/** Get argument rules */
const ArgumentRules& Func_MinBLTimeScaling::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "unscaledTree", Tree::getClassTypeSpec(), "The tree to be scaled.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa; has to match the tips of the tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "minBL", RealPos::getClassTypeSpec(), "Minimum branch length to use for time scaling.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_MinBLTimeScaling::getClassType(void)
{
    
    static std::string rev_type = "Func_MinBLTimeScaling";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_MinBLTimeScaling::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the primary Rev name for this function */
std::string Func_MinBLTimeScaling::getFunctionName( void ) const
{
    // create a name variable that is the same for all instances of this class
    std::string f_name = "fnMinBLTimeScaling";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_MinBLTimeScaling::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_MinBLTimeScaling::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TimeTree::getClassTypeSpec();
    return return_typeSpec;
}
