/**
 * @file
 * This file contains the declaration of the RevLanguage gtr function, which
 * is used to created deterministic variable associated with the gtr function.
 *
 * @brief Declaration and implementation of Func_gtr
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func__add.h 1406 2012-04-20 02:06:14Z hoehna $
 */


#ifndef Func_numUniqueInVector_H
#define Func_numUniqueInVector_H

#include "Natural.h"
#include "RlTypedFunction.h"

#include <string>
#include <cstdint>

namespace RevLanguage {
    
    template <typename valType>
    class Func_numUniqueInVector : public TypedFunction<Natural> {
        
    public:
        Func_numUniqueInVector( void );
        
        // Basic utility functions
        Func_numUniqueInVector*                         clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<std::int64_t>*               createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
       
    };
    
}

#include "Integer.h"
#include "ModelVector.h"
#include "NumUniqueInVector.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"

//using namespace RevLanguage;

/** default constructor */
template <typename valType>
RevLanguage::Func_numUniqueInVector<valType>::Func_numUniqueInVector( void ) : TypedFunction<Natural>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType>
RevLanguage::Func_numUniqueInVector<valType>* RevLanguage::Func_numUniqueInVector<valType>::clone( void ) const {
    
    return new RevLanguage::Func_numUniqueInVector<valType>( *this );
}

template <typename valType>
RevBayesCore::TypedFunction<std::int64_t>* RevLanguage::Func_numUniqueInVector<valType>::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<typename valType::valueType> >* vect = static_cast<const ModelVector<valType> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::NumUniqueInVector<typename valType::valueType>* f = new RevBayesCore::NumUniqueInVector<typename valType::valueType>( vect );
    
    return f;
}


/* Get argument rules */
template <typename valType>
const RevLanguage::ArgumentRules& RevLanguage::Func_numUniqueInVector<valType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "vector", ModelVector<valType>::getClassTypeSpec(), "The vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType>
const std::string& RevLanguage::Func_numUniqueInVector<valType>::getClassType(void)
{
    
    static std::string rev_type = "Func_numUniqueInVector";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_numUniqueInVector<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType>
std::string RevLanguage::Func_numUniqueInVector<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnNumUniqueInVector";
    
    return f_name;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_numUniqueInVector<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

#endif
