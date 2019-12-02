#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_simStartingTree.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RlClade.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "RlUtils.h"
#include "StartingTreeSimulator.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

#include <math.h>       /* log2 */

using namespace RevLanguage;

/** Default constructor */
Func_simStartingTree::Func_simStartingTree( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_simStartingTree* Func_simStartingTree::clone( void ) const
{
    
    return new Func_simStartingTree( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_simStartingTree::execute( void )
{
    
    const std::vector<RevBayesCore::Taxon> &taxa               = static_cast<const ModelVector<Taxon> &>( args[0].getVariable()->getRevObject() ).getValue();
    const RevBayesCore::RbVector<RevBayesCore::Clade>& constr  = static_cast<const ModelVector<Clade> &>( args[1].getVariable()->getRevObject() ).getValue();

    RevBayesCore::StartingTreeSimulator simulator;

    
    // the time tree object (topology + times)
    RevBayesCore::Tree *my_tree = simulator.simulateTree( taxa, constr );
    
    return new RevVariable( new TimeTree( my_tree ) );
}


/** Get argument rules */
const ArgumentRules& Func_simStartingTree::getArgumentRules( void ) const
{
    
    static ArgumentRules arg_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        arg_rules.push_back( new ArgumentRule( "taxa"  ,        ModelVector<Taxon>::getClassTypeSpec(), "The taxa used for initialization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        arg_rules.push_back( new ArgumentRule( "constraints",   ModelVector<Clade>::getClassTypeSpec(), "The topological constraints.",      ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<Clade>() ) );

        rules_set = true;
    }
    
    return arg_rules;
}


/** Get Rev type of object */
const std::string& Func_simStartingTree::getClassType(void)
{
    
    static std::string rev_type = "Func_simStartingTree";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_simStartingTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_simStartingTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "simStartingTree";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_simStartingTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_simStartingTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TimeTree::getClassTypeSpec();
    
    return return_typeSpec;
}
