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

#ifndef DistributionFunctionRv_H
#define DistributionFunctionRv_H

#include "RlTypedDistribution.h"
#include "Procedure.h"

#include <string>
#include <vector>

namespace RevLanguage {
    
    template <class valueType>
    class DistributionFunctionRv : public Procedure {
        
    public:
        DistributionFunctionRv(TypedDistribution<valueType> *d);                                                                //!< Object constructor
        DistributionFunctionRv(const DistributionFunctionRv& obj);                                                              //!< Copy constructor
        
        virtual                                 ~DistributionFunctionRv(void);
        
        // overloaded operators
        DistributionFunctionRv&                 operator=(const DistributionFunctionRv& c);
        
        // Basic utility functions
        DistributionFunctionRv*                 clone(void) const;                                                              //!< Clone the object
        static const std::string&               getClassType(void);                                                             //!< Get Rev type
        static const TypeSpec&                  getClassTypeSpec(void);                                                         //!< Get class type spec
        std::vector<std::string>                getFunctionNameAliases(void) const;                                             //!< Get the aliases of the name of the function in Rev
        std::string                             getFunctionName(void) const;                                                    //!< Get the primary name of the function in Rev
        const TypeSpec&                         getTypeSpec(void) const;                                                        //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                     execute(void);                                                                  //!< Execute the function. This is the function one has to overwrite for single return values.
        const ArgumentRules&                    getArgumentRules(void) const;                                                   //!< Get argument rules
        const TypeSpec&                         getReturnType(void) const;                                                      //!< Get type of return value
        
    protected:
        
        ArgumentRules                           argRules;                                                                       //!< Member rules converted to reference rules
        TypedDistribution<valueType>*           templateObject;                                                                 //!< The template object
        
    };
    
}

#include "ArgumentRule.h"
#include "DeterministicNode.h"
#include "DistributionFunctionRv.h"
#include "ModelVector.h"
#include "Natural.h"
#include "ProbabilityDensityFunction.h"
#include "Real.h"
#include "RevObject.h"
#include "TypeSpec.h"

#include <sstream>

/** Constructor */
template <class valueType>
RevLanguage::DistributionFunctionRv<valueType>::DistributionFunctionRv( TypedDistribution<valueType> *d ) : Procedure(),
    templateObject( d )
{
    
    argRules.push_back( new ArgumentRule("n", Natural::getClassTypeSpec(), "Number of random values to draw.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1)));
    const ArgumentRules &memberRules = templateObject->getParameterRules();
    for (std::vector<ArgumentRule*>::const_iterator it = memberRules.begin(); it != memberRules.end(); ++it)
    {
        argRules.push_back( (*it)->clone() );
    }
    
}


/** Constructor */
template <class valueType>
RevLanguage::DistributionFunctionRv<valueType>::DistributionFunctionRv(const DistributionFunctionRv& obj) : Procedure(obj),
    argRules( obj.argRules )
{
    
    templateObject = obj.templateObject->clone();
    
}


/** Destructor */
template <class valueType>
RevLanguage::DistributionFunctionRv<valueType>::~DistributionFunctionRv( void )
{
    
    delete templateObject;
    
}


template <class valueType>
RevLanguage::DistributionFunctionRv<valueType>& RevLanguage::DistributionFunctionRv<valueType>::operator=(const DistributionFunctionRv &c)
{
    
    if (this != &c)
    {
        Procedure::operator=(c);
        delete templateObject;
        
        templateObject = c.templateObject->clone();
        
        argRules = c.argRules;
    }
    
    return *this;
}


/** Clone the object */
template <class valueType>
RevLanguage::DistributionFunctionRv<valueType>* RevLanguage::DistributionFunctionRv<valueType>::clone(void) const
{
    
    return new DistributionFunctionRv<valueType>(*this);
}


/** Execute function: we reset our template object here and give out a copy */
template <class valueType>
RevLanguage::RevPtr<RevVariable> RevLanguage::DistributionFunctionRv<valueType>::execute( void )
{
    
    TypedDistribution<valueType>* copyObject = templateObject->clone();
    
    for ( size_t i = 1; i < args.size(); i++ )
    {
        
        if ( args[i].isConstant() )
        {
            copyObject->setConstParameter( args[i].getLabel(), RevPtr<const RevVariable>( (RevVariable*) args[i].getVariable() ) );
        }
        else
        {
            copyObject->setParameter( args[i].getLabel(), args[i].getReferenceVariable() );
        }
        
    }
    
    int n = (int)static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getValue();
    ModelVector<valueType> *values = new ModelVector<valueType>();
    for (int i = 0;  i < n; ++i)
    {
        valueType* value = copyObject->createRandomVariable();
        value->makeConstantValue();
        values->push_back( *value );
        delete value;
    }
    
    delete copyObject;
    
    return new RevVariable( values );
}


/** Get argument rules */
template <class valueType>
const RevLanguage::ArgumentRules& RevLanguage::DistributionFunctionRv<valueType>::getArgumentRules(void) const
{
    
    return argRules;
}


/** Get Rev type of object */
template <class valueType>
const std::string& RevLanguage::DistributionFunctionRv<valueType>::getClassType(void)
{
    
    static std::string rev_type = "DistributionFunctionRv";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
template <class valueType>
const RevLanguage::TypeSpec& RevLanguage::DistributionFunctionRv<valueType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the aliases for the function.
 * We simple return the aliases of the distribution.
 */
template <class valueType>
std::vector<std::string> RevLanguage::DistributionFunctionRv<valueType>::getFunctionNameAliases( void ) const
{
    
    std::vector<std::string> dist_aliases = ( templateObject != NULL ? templateObject->getDistributionFunctionAliases() : std::vector<std::string>() );
    std::vector<std::string> aliases;
    
    for (size_t i = 0; i < dist_aliases.size(); ++i)
    {
        std::string f_name = "r" + dist_aliases[i];
        aliases.push_back( f_name );
    }
    
    return aliases;
}


/**
 * Get the primary Rev name for this function.
 */
template <class valueType>
std::string RevLanguage::DistributionFunctionRv<valueType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "r" + templateObject->getDistributionFunctionName();
    
    return f_name;
}


/** Get type spec */
template <class valueType>
const RevLanguage::TypeSpec& RevLanguage::DistributionFunctionRv<valueType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
template <class valueType>
const RevLanguage::TypeSpec& RevLanguage::DistributionFunctionRv<valueType>::getReturnType(void) const
{
    
    return ModelVector< valueType >::getClassTypeSpec();
}


#endif
