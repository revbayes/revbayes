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
        
        std::vector<std::string>                                        getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                                     getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                                     getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                                     getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>                      getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                                        getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                                     getHelpTitle(void) const;                                   //!< Get the title of this help entry

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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
template <typename valType>
std::vector<std::string> RevLanguage::Func_ifelse<valType>::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
template <typename valType>
std::string RevLanguage::Func_ifelse<valType>::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "If the expression is true, then the function returns the first value, otherwise the second value.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
template <typename valType>
std::string RevLanguage::Func_ifelse<valType>::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details =  "The ifelse function is important when the value of a variable should deterministically change during an analysis depending on other variables. Standard if-else statements are not dynamically re-evaluated." ;
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
template <typename valType>
std::string RevLanguage::Func_ifelse<valType>::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "a <- 1\n";
    example += "b := ifelse( a == 1, 10, -10 )\n";
    example += "b\n";
    example += "\n";
    example += "a <- 2\n";
    example += "b\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
template <typename valType>
std::vector<RevBayesCore::RbHelpReference> RevLanguage::Func_ifelse<valType>::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
template <typename valType>
std::vector<std::string> RevLanguage::Func_ifelse<valType>::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
template <typename valType>
std::string RevLanguage::Func_ifelse<valType>::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "If-else statement as a function";
    
    return title;
}



/** Get type spec */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_ifelse<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif


