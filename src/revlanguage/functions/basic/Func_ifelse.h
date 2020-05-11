/**
 * @file
 * This file contains the declaration and implementation
 * of the templated Func_ifelse, which is used to create
 * value ifelses.
 *
 * @brief Declaration and implementation of Func_ifelse
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-12 10:25:58 +0200 (Tue, 12 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_ifelse.h 1626 2012-06-12 08:25:58Z hoehna $
 */

#ifndef Func_ifelse_H
#define Func_ifelse_H

#include "RlTypedFunction.h"

namespace RevLanguage {
    
    template <typename valType>
    class Func_ifelse : public TypedFunction<valType> {
        
    public:
        Func_ifelse();
        
        // Basic utility functions
        Func_ifelse*                                                    clone(void) const;                                          //!< Clone the object
        static const std::string&                                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                     getFunctionName(void) const;
        const TypeSpec&                                                 getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< typename valType::valueType>*      createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                            getArgumentRules(void) const;                               //!< Get argument rules
        
    protected:
        
    };
    
}


#include "ArgumentRule.h"
#include "RbUtil.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "IfElseFunction.h"



template <typename valType>
RevLanguage::Func_ifelse<valType>::Func_ifelse() : TypedFunction<valType>()
{
    
}

/* Clone object */
template <typename valType>
RevLanguage::Func_ifelse<valType>* RevLanguage::Func_ifelse<valType>::clone( void ) const
{
    
    return new Func_ifelse( *this );
}


template <typename valType>
RevBayesCore::TypedFunction< typename valType::valueType >* RevLanguage::Func_ifelse<valType>::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Boolean>       *c = static_cast<const RlBoolean &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<typename valType::valueType> *a = static_cast<const valType &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<typename valType::valueType> *b = static_cast<const valType &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::IfElseFunction<typename valType::valueType> *func = new RevBayesCore::IfElseFunction<typename valType::valueType>( c, a, b );
    
    return func;
}


/** Get argument rules */
template <typename valType>
const RevLanguage::ArgumentRules& RevLanguage::Func_ifelse<valType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "condition", RlBoolean::getClassTypeSpec(), "A variable representing the condition of the if-else statement.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "a"        , valType::getClassTypeSpec()  , "The value if the statement is true.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "b"        , valType::getClassTypeSpec()  , "The value if the statement is false.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Func_ifelse<valType>::getClassType(void)
{
    
    static std::string rev_type = "Func_ifelse<" + valType::getClassType() + ">";
    
	return rev_type;
}


/** Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_ifelse<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType>
std::string RevLanguage::Func_ifelse<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "ifelse";
    
    return f_name;
}



/** Get type spec */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_ifelse<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif


