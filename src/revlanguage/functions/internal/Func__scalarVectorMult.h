#ifndef Func__scalarVectorMult_H
#define Func__scalarVectorMult_H

#include "ModelVector.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the scalar-vector-Multiplication function.
     *
     * The RevLanguage wrapper of the scalar-vector-Multiplication function connects
     * the variables/parameters of the function and creates the internal ScalarVectorMultiplication object.
     * RevBayes > a <- [1,2,3]
     * RevBayes > b <- 2
     * RevBayes > a * b
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-31, version 1.0
     *
     */
    template <typename firstValType, typename secondValType, typename retType>
    class Func__scalarVectorMult :  public TypedFunction< retType > {
        
    public:
        Func__scalarVectorMult( void );
        
        // Basic utility functions
        Func__scalarVectorMult*                                         clone(void) const;                              //!< Clone the object
        static const std::string&                                       getClassType(void);                             //!< Get class name
        static const TypeSpec&                                          getClassTypeSpec(void);                         //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                        //!< Get the type spec of the instance
        bool                                                            isInternal(void) const { return true; }         //!< Is this an internal function?

        // Implementations of pure virtual functions of the base class(es)
        RevBayesCore::TypedFunction<typename retType::valueType>*       createFunction(void) const ;                    //!< Create a random variable from this distribution
        const ArgumentRules&                                            getArgumentRules(void) const;                   //!< Get argument rules
        
    private:
        
    };
    
}

#include "ScalarVectorMultiplication.h"
#include "RlDeterministicNode.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"

/** default constructor */
template <typename firstValType, typename secondValType, typename retType>
RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::Func__scalarVectorMult( void ) : TypedFunction<retType>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename firstValType, typename secondValType, typename retType>
RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>* RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::clone( void ) const
{
    
    return new Func__scalarVectorMult( *this );
}


template <typename firstValType, typename secondValType, typename retType>
RevBayesCore::TypedFunction<typename retType::valueType >* RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<typename firstValType::valueType>* firstArg = static_cast<const firstValType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<typename secondValType::valueType>* secondArg = static_cast<const secondValType &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::ScalarVectorMultiplication<typename firstValType::valueType, typename secondValType::elementType, typename retType::elementType> *func = new RevBayesCore::ScalarVectorMultiplication<typename firstValType::valueType, typename secondValType::elementType, typename retType::elementType>(firstArg, secondArg);
    
    return func;
}


/* Get argument rules */
template <typename firstValType, typename secondValType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "first" , firstValType::getClassTypeSpec() , "The left hand side variable." , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "second", secondValType::getClassTypeSpec(), "The right hand side variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename firstValType, typename secondValType, typename retType>
const std::string& RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::getClassType(void)
{
    static std::string rev_type = "Func__scalarVectorMult<" + firstValType::getClassType() + "," + secondValType::getClassType() + "," + retType::getClassType() + ">";
    
	return rev_type;
}

/* Get class type spec describing type of object */
template <typename firstValType, typename secondValType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( Func__scalarVectorMult<firstValType, secondValType, retType>::getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename firstValType, typename secondValType, typename retType>
std::string RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "mul";
    
    return f_name;
}


template <typename firstValType, typename secondValType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func__scalarVectorMult<firstValType, secondValType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

#endif
