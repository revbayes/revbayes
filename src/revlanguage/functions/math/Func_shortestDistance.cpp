//
//  Func_shortestDistance.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/15/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "Func_shortestDistance.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "ShortestDistanceFunction.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Natural.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_shortestDistance::Func_shortestDistance( void ) : TypedFunction<ModelVector<ModelVector<RealPos> > >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_shortestDistance* Func_shortestDistance::clone( void ) const
{
    
    return new Func_shortestDistance( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* Func_shortestDistance::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > >* adj = static_cast<const ModelVector<ModelVector<Natural> > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* dist = NULL;
    if ( this->args[1].getVariable()->getRevObject().isType( ModelVector<ModelVector<RealPos> >::getClassTypeSpec() ) )
    {
        dist = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    }
    else
    {
        size_t num_nodes = adj->getValue().size();
        dist = new RevBayesCore::ConstantNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >( "",
                                                                                                        new RevBayesCore::RbVector<RevBayesCore::RbVector<double> >(
                                                                                                                num_nodes,
                                                                                                                RevBayesCore::RbVector<double>(num_nodes, 1))
                                                                                                        );
    }
    
    RevBayesCore::ShortestDistanceFunction* f = new RevBayesCore::ShortestDistanceFunction( adj, dist );
    
    return f;

}


/* Get argument rules */
const ArgumentRules& Func_shortestDistance::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "adjacencies", ModelVector<ModelVector<Natural> >::getClassTypeSpec(), "The adjaceny matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "distances", ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "The distance matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        rules_set = true;
    }
    
    return argumentRules;
}

const std::string& Func_shortestDistance::getClassType(void)
{
    
    static std::string rev_type = "Func_shortestDistance";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_shortestDistance::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_shortestDistance::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnShortestDistance";
    
    return f_name;
}


const TypeSpec& Func_shortestDistance::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
