#ifndef Func__conversion_H
#define Func__conversion_H

#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    template <typename rlTypeFrom, typename rlTypeTo>
    class Func__conversion : public TypedFunction<rlTypeTo> {
        
    public:
        Func__conversion( void );
        
        // Basic utility functions
        Func__conversion*                                               clone(void) const;                              //!< Clone the object
        static const std::string&                                       getClassType(void);                             //!< Get class name
        static const TypeSpec&                                          getClassTypeSpec(void);                         //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                        //!< Get the type spec of the instance
        bool                                                            isInternal(void) const { return true; }         //!< Is this an internal function?

        // Implementations of pure virtual functions of the base class(es)
        RevBayesCore::TypedFunction<typename rlTypeTo::valueType>*      createFunction(void) const ;                    //!< Create a new internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                   //!< Get argument rules
        
    private:
        
    };
    
}

#include "TypeConversionFunction.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"

/** default constructor */
template <typename rlTypeFrom, typename rlTypeTo>
RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::Func__conversion( void ) : TypedFunction<rlTypeTo>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename rlTypeFrom, typename rlTypeTo>
RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>* RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::clone( void ) const
{
    
    return new Func__conversion( *this );
}


template <typename rlTypeFrom, typename rlTypeTo>
RevBayesCore::TypedFunction<typename rlTypeTo::valueType>* RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<typename rlTypeFrom::valueType>* arg = static_cast<const rlTypeFrom &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypeConversionFunction<typename rlTypeFrom::valueType, typename rlTypeTo::valueType> *func = new RevBayesCore::TypeConversionFunction<typename rlTypeFrom::valueType, typename rlTypeTo::valueType>(arg);
    
    return func;
}


/* Get argument rules */
template <typename rlTypeFrom, typename rlTypeTo>
const RevLanguage::ArgumentRules& RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "arg" , rlTypeFrom::getClassTypeSpec() , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename rlTypeFrom, typename rlTypeTo>
const std::string& RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::getClassType(void)
{
    static std::string rev_type = "Func__conversion<" + rlTypeFrom::getClassType() + "," + rlTypeTo::getClassType() + ">";
    
    return rev_type;
}

/* Get class type spec describing type of object */
template <typename rlTypeFrom, typename rlTypeTo>
const RevLanguage::TypeSpec& RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( Func__conversion<rlTypeFrom, rlTypeTo>::getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename rlTypeFrom, typename rlTypeTo>
std::string RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = rlTypeFrom::getClassType() + "2" + rlTypeTo::getClassType();
    
    return f_name;
}


template <typename rlTypeFrom, typename rlTypeTo>
const RevLanguage::TypeSpec& RevLanguage::Func__conversion<rlTypeFrom, rlTypeTo>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

#endif
