#ifndef Func_trunc_H
#define Func_trunc_H

#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the trunc function.
     *
     * The RevLanguage wrapper of the trunc function connects
     * the variables/parameters of the function and creates the internal TruncateFunction object.
     * Please read the TruncateFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    template <typename valType, typename retType>
    class Func_trunc : public TypedFunction<retType> {
        
    public:
        Func_trunc( void );
        
        // Basic utility functions
        Func_trunc*                                                             clone(void) const;                                          //!< Clone the object
        static const std::string&                                               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                         getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< typename retType::valueType>*              createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                    getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#include "TruncateFunction.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

/** default constructor */
template <typename valType, typename retType>
RevLanguage::Func_trunc<valType, retType>::Func_trunc( void ) : TypedFunction<retType>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType, typename retType>
RevLanguage::Func_trunc<valType, retType>* RevLanguage::Func_trunc<valType, retType>::clone( void ) const
{
    
    return new Func_trunc<valType, retType>( *this );
}


template <typename valType, typename retType>
RevBayesCore::TypedFunction< typename retType::valueType >* RevLanguage::Func_trunc<valType, retType>::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* arg = static_cast<const valType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TruncateFunction* f = new RevBayesCore::TruncateFunction( arg );
    
    return f;
}


/* Get argument rules */
template <typename valType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::Func_trunc<valType, retType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set ) 
    {
        
        argumentRules.push_back( new ArgumentRule( "x", valType::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType, typename retType>
const std::string& RevLanguage::Func_trunc<valType, retType>::getClassType(void)
{
    
    static std::string rev_type = "Func_trunc<" + valType::getClassType() + "," + retType::getClassType() + ">";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_trunc<valType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType, typename retType>
std::string RevLanguage::Func_trunc<valType, retType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "trunc";
    
    return f_name;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_trunc<valType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif
