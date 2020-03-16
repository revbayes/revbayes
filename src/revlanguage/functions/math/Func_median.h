#ifndef Func_median_H
#define Func_median_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the arithmetic median function.
     *
     * The RevLanguage wrapper of the median function connects
     * the variables/parameters of the function and creates the internal medianFunction object.
     *
     * @copybrief RevBayesCore::MedianFunction
     * @see RevBayesCore::MedianFunction for the internal object
     *
     */
    class Func_median : public TypedFunction<Real> {
        
    public:
        Func_median( void );
        
        // Basic utility functions
        Func_median*                                    clone(void) const;                                          //!< Clone the object
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
