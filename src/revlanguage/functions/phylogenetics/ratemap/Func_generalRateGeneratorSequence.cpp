//
//  Func_exp.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 8/7/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

//#include "BiogeographyRateGeneratorSequenceFunction.h"
#include "Func_generalRateGeneratorSequence.h"

#include <cstddef>

#include "GeneralRateGeneratorSequenceFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateGeneratorSequence.h"
#include "RlCharacterHistoryRateModifier.h"
#include "RlDeterministicNode.h"
#include "RlRateGeneratorSequence.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RateGenerator.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RlFunction.h"
#include "RlRateGenerator.h"
#include "TypeSpec.h"

namespace RevBayesCore { class CharacterHistoryRateModifier; }

using namespace RevLanguage;

/** default constructor */
Func_generalRateGeneratorSequence::Func_generalRateGeneratorSequence( void ) : TypedFunction<RateGeneratorSequence>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_generalRateGeneratorSequence* Func_generalRateGeneratorSequence::clone( void ) const {
    
    return new Func_generalRateGeneratorSequence( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RateGeneratorSequence>* Func_generalRateGeneratorSequence::createFunction() const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    size_t ns = rm->getValue().getNumberOfStates();
    unsigned long nc = static_cast<const Natural&>( this->args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::CharacterHistoryRateModifier> >* rm_vec = NULL;
    if ( this->args[2].getVariable()->getRevObject().isType( ModelVector<CharacterHistoryRateModifier>::getClassTypeSpec() ) )
    {
         rm_vec = static_cast<const ModelVector<CharacterHistoryRateModifier>&>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    }
    
    RevBayesCore::GeneralRateGeneratorSequenceFunction* f = new RevBayesCore::GeneralRateGeneratorSequenceFunction(ns, nc);
    f->setRateMatrix(rm);
    f->setRateModifiers(rm_vec);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_generalRateGeneratorSequence::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;
    
    if ( !rulesSet )
    {
        argumentRules.push_back( new ArgumentRule( "Q", RateGenerator::getClassTypeSpec(), "The per-character rate generator.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numChars", Natural::getClassTypeSpec(), "The number of characters.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rateModifiers", ModelVector<CharacterHistoryRateModifier>::getClassTypeSpec(), "The sequence-wide context-dependent rate modifiers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        
        rulesSet = true;
    }
    
    return argumentRules;
}


const std::string& Func_generalRateGeneratorSequence::getClassType(void)
{
    
    static std::string revType = "Func_generalRateGeneratorSequence";
    
	return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_generalRateGeneratorSequence::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return revTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_generalRateGeneratorSequence::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnRateGeneratorSequence";
    
    return f_name;
}


const TypeSpec& Func_generalRateGeneratorSequence::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/** Set a member variable */
void Func_generalRateGeneratorSequence::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "qSite" )
    {
        q = var;
    }
    else if ( name == "rateModifiers" )
    {
        rateModifiers = var;
    }
    else if ( name == "rfSite" )
    {
        root_frequencies = var;
    }
    else if ( name == "numChars" )
    {
        num_chars = var;
    }
    else
    {
        Function::setConstParameter(name, var);
    }
    
}
