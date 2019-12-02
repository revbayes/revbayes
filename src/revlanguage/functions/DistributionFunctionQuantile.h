/**
 * @file
 * This file contains the declaration of ConstructorFunction, which is used
 * for functions that construct objects using the RevLanguage.
 *
 * @brief Declaration of ConstructorFunction
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-01 14:40:37 +0200 (Fri, 01 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: ConstructorFunction.h 1602 2012-06-01 12:40:37Z hoehna $
 */

#ifndef DistributionFunctionQuantile_H
#define DistributionFunctionQuantile_H

#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"

#include <string>
#include <vector>

namespace RevLanguage {
    
    template <class valueType>
    class DistributionFunctionQuantile : public TypedFunction<valueType> {
        
    public:
        DistributionFunctionQuantile(TypedDistribution<valueType> *d);                                                                 //!< Object constructor
        DistributionFunctionQuantile(const DistributionFunctionQuantile& obj);                                                         //!< Copy constructor
        
        virtual                                         ~DistributionFunctionQuantile(void);
        
        // overloaded operators
        DistributionFunctionQuantile&                   operator=(const DistributionFunctionQuantile& c);
        
        // Basic utility functions
        DistributionFunctionQuantile*                   clone(void) const;                                                              //!< Clone the object
        static const std::string&                       getClassType(void);                                                             //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                         //!< Get class type spec
        std::vector<std::string>                        getFunctionNameAliases(void) const;                                             //!< Get the aliases of the name of the function in Rev
        std::string                                     getFunctionName(void) const;                                                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                                        //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                                     //!< Create a function object
        const ArgumentRules&                            getArgumentRules(void) const;                                                   //!< Get argument rules
        
    protected:
        
        ArgumentRules                                   argRules;                                                                               //!< Member rules converted to reference rules
        TypedDistribution<valueType>*                   templateObject;                                                                         //!< The template object
        
    };
    
}

#include "ArgumentRule.h"
#include "DeterministicNode.h"
#include "DistributionFunctionQuantile.h"
#include "Probability.h"
#include "ProbabilityDensityFunction.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "Real.h"
#include "TypeSpec.h"

#include <sstream>

/** Constructor */
template <class valueType>
RevLanguage::DistributionFunctionQuantile<valueType>::DistributionFunctionQuantile( TypedDistribution<valueType> *d ) : TypedFunction<Real>(),
templateObject( d )
{
    
    argRules.push_back( new ArgumentRule("x", valueType::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
    const ArgumentRules &memberRules = templateObject->getParameterRules();
    for (std::vector<ArgumentRule*>::const_iterator it = memberRules.begin(); it != memberRules.end(); ++it)
    {
        argRules.push_back( (*it)->clone() );
    }
    argRules.push_back( new ArgumentRule("log", RlBoolean::getClassTypeSpec(), ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
}


/** Constructor */
template <class valueType>
RevLanguage::DistributionFunctionQuantile<valueType>::DistributionFunctionQuantile(const DistributionFunctionQuantile& obj) : TypedFunction<Real>(obj),
argRules( obj.argRules )
{
    
    templateObject = obj.templateObject->clone();
    
}


/** Destructor */
template <class valueType>
RevLanguage::DistributionFunctionQuantile<valueType>::~DistributionFunctionQuantile( void )
{
    
    delete templateObject;
    
}


template <class valueType>
RevLanguage::DistributionFunctionQuantile<valueType>& RevLanguage::DistributionFunctionQuantile<valueType>::operator=(const DistributionFunctionQuantile &c) {
    
    if (this != &c) {
        Function::operator=(c);
        
        delete templateObject;
        templateObject = c.templateObject->clone();
        
        argRules = c.argRules;
    }
    
    return *this;
}


/** Clone the object */
template <class valueType>
RevLanguage::DistributionFunctionQuantile<valueType>* RevLanguage::DistributionFunctionQuantile<valueType>::clone(void) const
{
    
    return new DistributionFunctionQuantile<valueType>(*this);
}


/** Execute function: we reset our template object here and give out a copy */
template <class valueType>
RevBayesCore::TypedFunction<double>* RevLanguage::DistributionFunctionQuantile<valueType>::createFunction( void ) const
{
    
    TypedDistribution<valueType>* copyObject = templateObject->clone();
    
    for ( size_t i = 1; i < (args.size()-1); i++ )
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
    
    
    
    RevBayesCore::TypedDagNode<typename valueType::valueType>* arg = static_cast<const valueType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::ProbabilityDensityFunction<typename valueType::valueType>* f = new RevBayesCore::ProbabilityDensityFunction<typename valueType::valueType>( arg, copyObject->createDistribution() );
    
    return f;
}


/** Get argument rules */
template <class valueType>
const RevLanguage::ArgumentRules& RevLanguage::DistributionFunctionQuantile<valueType>::getArgumentRules(void) const {
    
    return argRules;
}


/** Get Rev type of object */
template <class valueType>
const std::string& RevLanguage::DistributionFunctionQuantile<valueType>::getClassType(void)
{
    
    static std::string rev_type = "DistributionFunctionQuantile";
    
    return rev_type;
}


/** Get class type spec describing type of object */
template <class valueType>
const RevLanguage::TypeSpec& RevLanguage::DistributionFunctionQuantile<valueType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the aliases for the function.
 * We simple return the aliases of the distribution.
 */
template <class valueType>
std::vector<std::string> RevLanguage::DistributionFunctionQuantile<valueType>::getFunctionNameAliases( void ) const
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


/**
 * Get the primary Rev name for this function.
 */
template <class valueType>
std::string RevLanguage::DistributionFunctionQuantile<valueType>::getFunctionName( void ) const
{
    // create a name variable that is NOT the same for all instance of this class
    std::string f_name = "q" + templateObject->getDistributionFunctionName();
    
    return f_name;
}


/** Get type spec */
template <class valueType>
const RevLanguage::TypeSpec& RevLanguage::DistributionFunctionQuantile<valueType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif
