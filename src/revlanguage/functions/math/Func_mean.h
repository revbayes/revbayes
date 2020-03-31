#ifndef Func_mean_H
#define Func_mean_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the arithmetic mean function.
     *
     * The RevLanguage wrapper of the mean function connects
     * the variables/parameters of the function and creates the internal MeanFunction object.
     *
     * @copybrief RevBayesCore::MeanFunction
     * @see RevBayesCore::MeanFunction for the internal object
     *
     *
     *
     */
    class Func_mean : public TypedFunction<Real> {
        
    public:
        Func_mean( void );
        
        // Basic utility functions
        Func_mean*                                      clone(void) const;                                          //!< Clone the object
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
