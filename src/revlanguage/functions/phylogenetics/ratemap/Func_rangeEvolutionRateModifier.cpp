//
//  Func_rangeEvolutionRateModifier.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "Func_rangeEvolutionRateModifier.h"

#include <cstddef>
#include <ostream>

#include "RangeEvolutionRateModifierFunction.h"
#include "CharacterHistoryRateModifier.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbException.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_rangeEvolutionRateModifier::Func_rangeEvolutionRateModifier( void ) : TypedFunction< CharacterHistoryRateModifier >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_rangeEvolutionRateModifier* Func_rangeEvolutionRateModifier::clone( void ) const {
    
    return new Func_rangeEvolutionRateModifier( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >* Func_rangeEvolutionRateModifier::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* gf = static_cast<const Real&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* lf = static_cast<const Real&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* c = NULL;
    if ( this->args[2].getVariable()->getRevObject().isType( ModelVector<ModelVector<Real> >::getClassTypeSpec() ) )
    {
        c = static_cast<const ModelVector<ModelVector<Real> > &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    }
    else
    {
        throw RbException("Failed to type cast `matrix` argument!");
    }
    size_t nc = static_cast<const Natural&>( this->args[3].getVariable()->getRevObject() ).getValue();
    
    
    // validate input
    if (c != NULL) {
        size_t d1 = c->getValue().size();
        size_t d2 = c->getValue()[0].size();
        if (d1 != nc || d2 != nc)
        {
            std::stringstream ss;
            ss << "Context matrix dimensions (" << d1 << " x " << d2 << ") do not match number of characters (" << nc << ")";
            throw RbException(ss.str());
        }
    }
    
    
    RevBayesCore::RangeEvolutionRateModifierFunction* f = new RevBayesCore::RangeEvolutionRateModifierFunction(gf, lf, c, nc);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_rangeEvolutionRateModifier::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;
    
    if ( !rulesSet )
    {
        
        argumentRules.push_back( new ArgumentRule( "gainFactor", Real::getClassTypeSpec(), "Multiplicative factor (r' = r * e^{ n_1 * f }) for characters in context set", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "lossFactor", Real::getClassTypeSpec(), "Multiplicative factor (r' = r * e^{ n_0 * f }) for characters in context set", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "matrix",  ModelVector<ModelVector<Real> >::getClassTypeSpec(), "Weighted character adjacency matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "numChars", Natural::getClassTypeSpec(), "Number of characters", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        
        rulesSet = true;
    }
    
    return argumentRules;
}


const std::string& Func_rangeEvolutionRateModifier::getClassType(void)
{
    
    static std::string revType = "Func_rangeEvolutionRateModifier";
    
    return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_rangeEvolutionRateModifier::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_rangeEvolutionRateModifier::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnRangeEvolutionRateModifier";
    
    return f_name;
}


const TypeSpec& Func_rangeEvolutionRateModifier::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}