/**
 * @file
 * This file contains the declaration and implementation
 * of the templated Func_vector, which is used to create
 * value vectors.
 *
 * @brief Declaration and implementation of Func_vector
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-12 10:25:58 +0200 (Tue, 12 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_vector.h 1626 2012-06-12 08:25:58Z hoehna $
 */

#ifndef Func__eq_H
#define Func__eq_H

#include "RlBoolean.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    template <typename leftValType, typename rightValType>
    class Func__eq : public TypedFunction<RlBoolean> {
        
    public:
        Func__eq();
        
        // Basic utility functions
        Func__eq*                                               clone(void) const;                                          //!< Clone the object
        static const std::string&                               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                         getTypeSpec(void) const;                                    //!< Get language type of the object
        bool                                                    isInternal(void) const { return true; }                     //!< Is this an internal function?

        // Regular functions
        RevBayesCore::TypedFunction<RevBayesCore::Boolean>*     createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                    getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#include "ArgumentRule.h"
#include "EquationFunction.h"
#include "RbUtil.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"



template <typename leftValType, typename rightValType>
RevLanguage::Func__eq<leftValType,rightValType>::Func__eq() : TypedFunction<RlBoolean>()
{
    
}

/* Clone object */
template <typename leftValType, typename rightValType>
RevLanguage::Func__eq<leftValType,rightValType>* RevLanguage::Func__eq<leftValType,rightValType>::clone( void ) const
{
    
    return new Func__eq<leftValType,rightValType>( *this );
}


/** Execute function: We rely on getValue and overloaded push_back to provide functionality */
template <typename leftValType, typename rightValType>
RevBayesCore::TypedFunction<RevBayesCore::Boolean>* RevLanguage::Func__eq<leftValType,rightValType>::createFunction( void ) const
{
        
    const RevBayesCore::TypedDagNode<typename leftValType::valueType>* leftVal = static_cast<const leftValType &>( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<typename rightValType::valueType>* rightVal = static_cast<const rightValType &>( args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::EquationFunction<typename leftValType::valueType, typename rightValType::valueType> *func = new RevBayesCore::EquationFunction<typename leftValType::valueType, typename rightValType::valueType>( leftVal, rightVal );
    
    return func;
}


/** Get argument rules */
template <typename leftValType, typename rightValType>
const RevLanguage::ArgumentRules& RevLanguage::Func__eq<leftValType,rightValType>::getArgumentRules( void ) const {
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "", leftValType::getClassTypeSpec() , "The left hand side of the expression.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "", rightValType::getClassTypeSpec(), "The right hand side of the expression.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
template <typename leftValType, typename rightValType>
const std::string& RevLanguage::Func__eq<leftValType,rightValType>::getClassType(void)
{
    
    static std::string rev_type = "Func__eq<" + leftValType::getClassType() + "," + rightValType::getClassType() + ">";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
template <typename leftValType, typename rightValType>
const RevLanguage::TypeSpec& RevLanguage::Func__eq<leftValType,rightValType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
template <typename leftValType, typename rightValType>
std::string RevLanguage::Func__eq<leftValType, rightValType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "eq";
    
    return f_name;
}


/** Get type spec */
template <typename leftValType, typename rightValType>
const RevLanguage::TypeSpec& RevLanguage::Func__eq<leftValType,rightValType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif


