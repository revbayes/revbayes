#ifndef Func_inverse_h
#define Func_inverse_h

#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the inverse function (inverse()).
     *
     * When passed 
     *
     *
     * @copyright Copyright 2024-
     * @author Martin R. Smith
     * @since 2024-07-11, version 1.2.5
     *
     */
    class Func_inverse :  public TypedFunction<Real> {
        
    public:
        // Basic utility functions
        Func_inverse*                                   clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif /* Func_inverse_h */
