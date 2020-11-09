#ifndef Func_min_H
#define Func_min_H

#include "RlTypedFunction.h"
#include "Real.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the minimum value function.
     *
     * The RevLanguage wrapper of the sum function connects
     * the variables/parameters of the function and creates the internal MinFunction object.
     *
     * @copybrief RevBayesCore::MinFunction
     * @see RevBayesCore::MinFunction for the internal object
     */
    class Func_min :  public TypedFunction<Real> {
        
    public:
        Func_min( void );
        
        // Basic utility functions
        Func_min*                                       clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
