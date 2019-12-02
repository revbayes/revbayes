//
//  Func_EpochCladoProbs.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 10/24/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "Func_EpochCladoProbs.h"
#include "ConstantNode.h"
#include "EpochCladogeneticStateFunction.h"
#include "MatrixReal.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "Real.h"
#include "RealPos.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlDeterministicNode.h"
#include "RlSimplex.h"
#include "RlString.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_EpochCladoProbs::Func_EpochCladoProbs( void ) : TypedFunction<CladogeneticProbabilityMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_EpochCladoProbs* Func_EpochCladoProbs::clone( void ) const {
    
    return new Func_EpochCladoProbs( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticProbabilityMatrix >* Func_EpochCladoProbs::createFunction( void ) const
{
    
    // supplied arguments
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::CladogeneticProbabilityMatrix> >* cp = static_cast<const ModelVector<CladogeneticProbabilityMatrix> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* et = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    if (et->getValue().size() != cp->getValue().size())
    {
        throw RbException("cladogeneticProbabilities and times must be of equal size.");
    }
    
    // create P matrix
    RevBayesCore::EpochCladogeneticStateFunction* f = new RevBayesCore::EpochCladogeneticStateFunction( et, cp, 0, 0 );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_EpochCladoProbs::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "cladogeneticProbabilities",
                                                  ModelVector<CladogeneticProbabilityMatrix>::getClassTypeSpec(),
                                                  "The cladogenetic event probability functions.",
                                                  ArgumentRule::BY_CONSTANT_REFERENCE,
                                                  ArgumentRule::ANY ));
        
        argumentRules.push_back( new ArgumentRule( "times", ModelVector<RealPos>::getClassTypeSpec(), "The vector of epoch end ages.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}



const std::string& Func_EpochCladoProbs::getClassType(void)
{
    
    static std::string rev_type = "Func_EpochCladoProbs";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_EpochCladoProbs::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_EpochCladoProbs::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnEpochCladoProbs";
    
    return f_name;
}

const TypeSpec& Func_EpochCladoProbs::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
