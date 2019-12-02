#ifndef Func_hyperbolicSine_h
#define Func_hyperbolicSine_h

#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the hyperbolic Sine function (sinh()).
     *
     * The RevLanguage wrapper of the hyperbolic function function connects
     * the variables/parameters of the function and creates the internal HyperbolicSineFunction object.
     * Please read the HyperbolicSineFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Wade Dismukes)
     * @since 2016-09-26, version 1.0
     *
     */
    class Func_hyperbolicSine :  public TypedFunction<Real> {
        
    public:
        Func_hyperbolicSine( void );
        
        // Basic utility functions
        Func_hyperbolicSine*                         clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif /* Func_hyperbolicSine_h */
