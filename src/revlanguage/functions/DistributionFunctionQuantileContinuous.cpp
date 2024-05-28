//
//  DistributionFunctionQuantileContinuous.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 3/8/13.
//  Copyright 2013 __MyCompanyName__. All rights reserved.
//

#include "DistributionFunctionQuantileContinuous.h"

#include <cstddef>

#include "ArgumentRule.h"
#include "DeterministicNode.h"
#include "Probability.h"
#include "QuantileFunction.h"
#include "Real.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "RlFunction.h"

namespace RevBayesCore { class ContinuousDistribution; }

using namespace RevLanguage;

/** Constructor */
DistributionFunctionQuantileContinuous::DistributionFunctionQuantileContinuous( ContinuousDistribution *d ) : TypedFunction<Real>(),
    templateObject( d )
{
    
    argRules.push_back( new ArgumentRule("p", Probability::getClassTypeSpec(), "The probability (i.e., quantile) of the distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    const ArgumentRules &memberRules = templateObject->getParameterRules();
    for (std::vector<ArgumentRule*>::const_iterator it = memberRules.begin(); it != memberRules.end(); ++it)
    {
        argRules.push_back( (*it)->clone() );
    }
    
}


/** Constructor */
DistributionFunctionQuantileContinuous::DistributionFunctionQuantileContinuous(const DistributionFunctionQuantileContinuous& obj) : TypedFunction<Real>(obj),
    argRules( obj.argRules )
{
    
    if ( obj.templateObject != NULL )
    {
        templateObject = obj.templateObject->clone();
    }
    else
    {
        templateObject = NULL;
    }
    
}



DistributionFunctionQuantileContinuous::~DistributionFunctionQuantileContinuous( void )
{
    delete templateObject;
}

DistributionFunctionQuantileContinuous& DistributionFunctionQuantileContinuous::operator=(const DistributionFunctionQuantileContinuous &c) {
    
    if (this != &c)
    {
        TypedFunction<Real>::operator=(c);
        
        delete templateObject;
        if ( c.templateObject != NULL )
        {
            templateObject = c.templateObject->clone();
        }
        else
        {
            templateObject = NULL;
        }
        
        
        argRules = c.argRules;
    }
    
    return *this;
}


/** Clone the object */
DistributionFunctionQuantileContinuous* DistributionFunctionQuantileContinuous::clone(void) const
{
    
    return new DistributionFunctionQuantileContinuous(*this);
}


RevBayesCore::TypedFunction<double>* DistributionFunctionQuantileContinuous::createFunction( void ) const
{
    
    RevBayesCore::TypedFunction<double>* f = NULL;
    RevBayesCore::TypedDagNode<double>* arg = static_cast<const Probability &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    ContinuousDistribution* copyObject = templateObject->clone();
    
    for ( size_t i = 1; i < args.size(); i++ )
    {
    
        if ( args[i].isConstant() )
        {
            copyObject->setConstParameter( args[i].getLabel(), RevPtr<const RevVariable>( (const RevVariable*) args[i].getVariable() ) );
        }
        else
        {
            copyObject->setParameter( args[i].getLabel(), args[i].getReferenceVariable() );
        }
    }
        
    RevBayesCore::ContinuousDistribution *d = copyObject->createDistribution();
    f = new RevBayesCore::QuantileFunction( arg, d );
    
    
    // return the value
    return f;
}


/** Get argument rules */
const ArgumentRules& DistributionFunctionQuantileContinuous::getArgumentRules(void) const
{
    
    return argRules;
}


/** Get Rev type of object */
const std::string& DistributionFunctionQuantileContinuous::getClassType(void)
{
    
    static std::string rev_type = "DistributionFunctionQuantileContinuous";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& DistributionFunctionQuantileContinuous::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string DistributionFunctionQuantileContinuous::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "q" + templateObject->getDistributionFunctionName();
    
    return f_name;
}


/**
 * Get the aliases for the function.
 * We simple return the aliases of the distribution.
 */
std::vector<std::string> DistributionFunctionQuantileContinuous::getFunctionNameAliases( void ) const
{
    
    std::vector<std::string> dist_aliases = ( templateObject != NULL ? templateObject->getDistributionFunctionAliases() : std::vector<std::string>() );
    std::vector<std::string> aliases;
    
    for (size_t i = 0; i < dist_aliases.size(); ++i)
    {
        std::string f_name = "q" + dist_aliases[i];
        aliases.push_back( f_name );
    }
    
    return aliases;
}


/** Get type spec */
const TypeSpec& DistributionFunctionQuantileContinuous::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
