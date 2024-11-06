#ifndef Func__vectorIndexOperator_H
#define Func__vectorIndexOperator_H

#include "RlTypedFunction.h"

#include <string>
#include <cstdint>

namespace RevLanguage {
    
    template <typename valType>
    class Func__vectorIndexOperator : public TypedFunction<valType> {
        
    public:
        Func__vectorIndexOperator( void );
        
        // Basic utility functions
        Func__vectorIndexOperator*                                      clone(void) const;                              //!< Clone the object
        static const std::string&                                       getClassType(void);                             //!< Get class name
        static const TypeSpec&                                          getClassTypeSpec(void);                         //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                        //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<typename valType::valueType>*       createFunction(void) const ;                    //!< Create a new internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                   //!< Get argument rules
        
    private:
        
    };
    
}

#include "RlDeterministicNode.h"
#include "ModelVector.h"
#include "VectorIndexOperator.h"
#include "TypedDagNode.h"

/** default constructor */
template <typename valType>
RevLanguage::Func__vectorIndexOperator<valType>::Func__vectorIndexOperator( void ) : TypedFunction<valType>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename valType>
RevLanguage::Func__vectorIndexOperator<valType>* RevLanguage::Func__vectorIndexOperator<valType>::clone( void ) const
{
    
    return new Func__vectorIndexOperator( *this );
}


template <typename valType>
RevBayesCore::TypedFunction< typename valType::valueType >* RevLanguage::Func__vectorIndexOperator<valType>::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<typename valType::valueType> >* v = static_cast<const ModelVector<valType> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* index = static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::VectorIndexOperator<typename valType::valueType> *func = new RevBayesCore::VectorIndexOperator<typename valType::valueType>(v, index);
    
    return func;
}


/* Get argument rules */
template <typename valType>
const RevLanguage::ArgumentRules& RevLanguage::Func__vectorIndexOperator<valType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "v"    , ModelVector<valType>::getClassTypeSpec(), "The vector.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "index", Natural::getClassTypeSpec()             , "The index.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType>
const std::string& RevLanguage::Func__vectorIndexOperator<valType>::getClassType(void)
{
    
    static std::string revClassType = "Func__vectorIndexOperator<" + valType::getClassType() + ">";
    
	return revClassType;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func__vectorIndexOperator<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return revClassTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename valType>
std::string RevLanguage::Func__vectorIndexOperator<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "[]";
    
    return f_name;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func__vectorIndexOperator<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

#endif
