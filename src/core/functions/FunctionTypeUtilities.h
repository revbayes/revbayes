#ifndef FunctionTypeUtilities_H
#define FunctionTypeUtilities_H

#include "Function.h"
#include "RbException.h"
#include "TypedFunction.h"

namespace RevBayesCore {

    /** Recover a function's known result type, checking the invariant in debug builds. */
    template <typename valueType>
    TypedFunction<valueType>* assumeFunctionReturns(Function *function)
    {
#ifndef NDEBUG
        TypedFunction<valueType> *typedFunction = dynamic_cast<TypedFunction<valueType>*>( function );
        if ( typedFunction == NULL )
        {
            throw RbException("A function does not return the expected value type.");
        }
        return typedFunction;
#else
        return static_cast<TypedFunction<valueType>*>( function );
#endif
    }

    /** Recover a function's known result type, checking the invariant in debug builds. */
    template <typename valueType>
    const TypedFunction<valueType>* assumeFunctionReturns(const Function *function)
    {
#ifndef NDEBUG
        const TypedFunction<valueType> *typedFunction = dynamic_cast<const TypedFunction<valueType>*>( function );
        if ( typedFunction == NULL )
        {
            throw RbException("A function does not return the expected value type.");
        }
        return typedFunction;
#else
        return static_cast<const TypedFunction<valueType>*>( function );
#endif
    }

}

#endif
