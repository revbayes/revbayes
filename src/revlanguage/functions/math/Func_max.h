#ifndef Func_max_H
#define Func_max_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the max value function.
     *
     * The RevLanguage wrapper of the sum function connects
     * the variables/parameters of the function and creates the internal MaxFunction object.
     *
     * @copybrief RevBayesCore::MaxFunction
     * @see RevBayesCore::MaxFunction for the internal object
     */
    class Func_max :  public TypedFunction<Real> {
        
    public:
        Func_max( void );
        
        // Basic utility functions
        Func_max*                                       clone(void) const;                                          //!< Clone the object
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
