/**
 * @file
 * This file contains the declaration and implementation
 * of the templated Func__uminus, which is used to take the negative value of a variable.
 *
 * @brief Declaration and implementation of Func__add
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func__add.h 1406 2012-04-20 02:06:14Z hoehna $
 */

#ifndef Func_uminus_H
#define Func_uminus_H

#include "TypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {

template <typename firstValType, typename retType>
class Func__uminus : public TypedFunction<retType> {
    
public:
    Func__uminus( void );
    
    // Basic utility functions
    Func__uminus*                                   clone(void) const;                                          //!< Clone the object
    static const std::string&                       getClassType(void);                                         //!< Get Rev type
    static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
    std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
    const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
    bool                                            isInternal(void) const { return true; }                     //!< Is this an internal function?

    // Function functions you have to override
    RevBayesCore::TypedFunction<typename retType::valueType>*     createFunction(void) const;                   //!< Create a function object
    const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
    
};
    
}

#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "UnaryMinus.h"

/** default constructor */
template <typename firstValType, typename retType>
RevLanguage::Func__uminus<firstValType, retType>::Func__uminus( void ) : TypedFunction<retType>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename firstValType, typename retType>
RevLanguage::Func__uminus<firstValType, retType>* RevLanguage::Func__uminus<firstValType, retType>::clone( void ) const {
    
    return new Func__uminus( *this );
}


template <typename firstValType, typename retType>
RevBayesCore::TypedFunction<typename retType::valueType>* RevLanguage::Func__uminus<firstValType, retType>::createFunction(void) const
{
    
    RevBayesCore::TypedDagNode<typename firstValType::valueType>* firstArg = static_cast<const firstValType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::UnaryMinus<typename firstValType::valueType> *func = new RevBayesCore::UnaryMinus<typename firstValType::valueType>(firstArg);

    return func;
}


/* Get argument rules */
template <typename firstValType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::Func__uminus<firstValType, retType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set ) {
        
        argumentRules.push_back( new ArgumentRule( "value", firstValType::getClassTypeSpec(), "The variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename firstValType, typename retType>
const std::string& RevLanguage::Func__uminus<firstValType, retType>::getClassType(void)
{
    
    static std::string rev_type = "Func__uminus<" + firstValType::getClassType() + "," + retType::getClassType() + ">";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
template <typename firstValType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func__uminus<firstValType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
template <typename firstValType, typename retType>
std::string RevLanguage::Func__uminus<firstValType, retType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "uminus";
    
    return f_name;
}


template <typename firstValType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func__uminus<firstValType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

#endif
