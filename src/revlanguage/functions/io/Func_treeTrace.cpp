//
//  Func_TreeTrace.cpp
//  RevBayes_development_branch
//
//  Created by will freyman on 12/15/15.
//  Copyright (c) 2015 will freyman. All rights reserved.
//



#include "Func_treeTrace.h"
#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RbException.h"
#include "RlTree.h"
#include "RlTraceTree.h"
#include "TypedDagNode.h"

#include <map>
#include <set>
#include <sstream>


using namespace RevLanguage;


Func_treeTrace* Func_treeTrace::clone( void ) const
{
    
    return new Func_treeTrace( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_treeTrace::execute( void )
{
    RevBayesCore::TraceTree t(false);
    
    RevObject& ro = args[0].getVariable()->getRevObject();

    const ModelVector<RevLanguage::Tree>& trees = static_cast<const ModelVector<RevLanguage::Tree>&>( ro );

    if ( ro.isType( ModelVector<TimeTree>::getClassTypeSpec() ) )
    {
        t = RevBayesCore::TraceTree( true );
    }

    for (size_t i = 0; i < trees.size(); ++i)
    {
        t.addObject( new RevBayesCore::Tree( trees[i] ) );
    }

    int burnin = 0;

    RevObject& b = args[1].getVariable()->getRevObject();
    if ( b.isType( Integer::getClassTypeSpec() ) )
    {
        burnin = (int)static_cast<const Integer &>(b).getValue();
    }
    else
    {
        double burninFrac = static_cast<const Probability &>(b).getValue();
        burnin = int( floor( trees.size()*burninFrac ) );
    }

    t.setBurnin(burnin);

    return new RevVariable( new TraceTree( t ) );
}


/** Get argument rules */
const ArgumentRules& Func_treeTrace::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        std::vector<TypeSpec> treeTypes;
        treeTypes.push_back( ModelVector<TimeTree>::getClassTypeSpec() );
        treeTypes.push_back( ModelVector<BranchLengthTree>::getClassTypeSpec() );
        treeTypes.push_back( ModelVector<Tree>::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "trees", treeTypes, "Vector of trees.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> burninTypes;
        burninTypes.push_back( Probability::getClassTypeSpec() );
        burninTypes.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burninTypes     , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_treeTrace::getClassType(void)
{
    
    static std::string rev_type = "Func_treeTrace";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_treeTrace::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_treeTrace::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Will Freyman" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_treeTrace::getHelpDescription(void) const
{
    std::string description = "";
    description += "Creates a tree trace object from a vector of trees.";
    
    return description;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_treeTrace::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# Read in a vector of trees\n";
    example += "trees = readTrees(\"trees.nex\")\n";
    example += "\n";
    example += "# Create a tree trace\n";
    example += "tree_trace = treeTrace(trees, burnin=0.25)\n";
    example += "\n";
    example += "# Create a distribution of trees from the tree trace\n";
    example += "tree ~ dnEmpiricalTree(tree_trace)\n";
    example += "\n";
    example += "# Add an MCMC move\n";
    example += "moves[1] = mvEmpiricalTree(tree)\n";
    
    return example;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_treeTrace::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "mvEmpiricalTree" );
    see_also.push_back( "treeTrace" );
    see_also.push_back( "readTreeTrace" );
    see_also.push_back( "readTrees" );
    
    
    return see_also;
}




/**
 * Get the primary Rev name for this function.
 */
std::string Func_treeTrace::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "treeTrace";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_treeTrace::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_treeTrace::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TraceTree::getClassTypeSpec();
    return return_typeSpec;
}
