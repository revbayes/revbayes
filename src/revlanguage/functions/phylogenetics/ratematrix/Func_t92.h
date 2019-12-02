#ifndef Func_t92_H
#define Func_t92_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Tamura (1992) rate matrix function.
     *
     * The RevLanguage wrapper of the Tamura two parameter rate matrix connects
     * the variables/parameters of the function and creates the internal T92RateMatrixFunction object.
     * Please read the T92RateMatrixFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot and Sebastian Hoehna)
     * @since 2014-08-14, version 1.0
     *
     */
    class Func_t92 : public TypedFunction<RateMatrix> {
        
    public:
        Func_t92( void );
        
        // Basic utility functions
        Func_t92*                                                           clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >*         createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
