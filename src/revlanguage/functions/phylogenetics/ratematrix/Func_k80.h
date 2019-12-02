#ifndef Func_k80_H
#define Func_k80_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Kimura (1980) rate matrix function.
     *
     * The RevLanguage wrapper of the Kimura one parameter rate matrix connects
     * the variables/parameters of the function and creates the internal K80RateMatrixFunction object.
     * Please read the K80RateMatrixFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-14, version 1.0
     *
     */
    class Func_k80 : public TypedFunction<RateMatrix> {
        
    public:
        Func_k80( void );
        
        // Basic utility functions
        Func_k80*                                                           clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::vector<std::string>                                            getFunctionNameAliases( void ) const;                       //!< Get aliases of the function name
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >*         createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
