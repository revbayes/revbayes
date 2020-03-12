#ifndef Func_exp_H
#define Func_exp_H

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the exp function.
     *
     * The RevLanguage wrapper of the exp function connects
     * the variables/parameters of the function and creates the internal ExponentialFunction object.
     *
     * @copybrief RevBayesCore::ExponentialFunction
     * @see RevBayesCore::ExponentialFunction for the internal object
     */
    class Func_exp : public TypedFunction<RealPos> {
        
    public:
        Func_exp( void );
        
        // Basic utility functions
        Func_exp*                                       clone(void) const;                                          //!< Clone the object
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
