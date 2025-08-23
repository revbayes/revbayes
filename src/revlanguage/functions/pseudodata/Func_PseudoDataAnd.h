#ifndef Func_PseudoDataAnd_H
#define Func_PseudoDataAnd_H

#include <string>
#include <iosfwd>
#include <vector>

#include "Real.h"
#include "RlPseudoData.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;

    /**
     * The RevLanguage wrapper of And pseudo-data.
     *
     * @copyright Copyright 2025-
     * @author Benjamin D. Redelings
     * @since 2025-05-30, version 1.0
     *
     */
    template <typename T>
    class Func_PseudoDataAnd: public TypedFunction<PseudoData<T>> {

    public:
        Func_PseudoDataAnd( void ): TypedFunction<PseudoData<T>>( ) { }

        // Basic utility functions
        Func_PseudoDataAnd<T>* clone( void ) const
        {
            return new Func_PseudoDataAnd<T>( *this );
        }

        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::PseudoData<typename T::valueType>>*     createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#include "RlDeterministicNode.h"
#include "RlPseudoData.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlBoolean.h"
#include "Integer.h"
#include "Simplex.h"
#include "TypeSpec.h"

using std::vector;
using std::unique_ptr;

namespace Core = RevBayesCore;

template <typename T>
Core::PseudoData<T>* PseudoDataAndFunc(const Core::PseudoData<T>& left, const Core::PseudoData<T>& right) 
{
    typename Core::PseudoData<T>::func_t and_func = [=](const T& x)
        {
            return left(x) + right(x);
        };
    return new Core::PseudoData<T>(and_func);
}

using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */


template <typename T>
Core::TypedFunction< Core::PseudoData<typename T::valueType> >* Func_PseudoDataAnd<T>::createFunction( void ) const
{
    auto left  = dynamic_cast<const PseudoData<T> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    auto right = dynamic_cast<const PseudoData<T> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoData<typename T::valueType> >( PseudoDataAndFunc<typename T::valueType>, left, right );
}


/* Get argument rules */
template <typename T>
const ArgumentRules& Func_PseudoDataAnd<T>::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "left", PseudoData<T>::getClassTypeSpec(), "The left argument.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "right", PseudoData<T>::getClassTypeSpec(), "The right argument.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


template <typename T>
const std::string& Func_PseudoDataAnd<T>::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataAnd";

    return rev_type;
}

/* Get class type spec describing type of object */
template <typename T>
const TypeSpec& Func_PseudoDataAnd<T>::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename T>
std::string Func_PseudoDataAnd<T>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "_and";

    return f_name;
}


template <typename T>
const TypeSpec& Func_PseudoDataAnd<T>::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
#endif


