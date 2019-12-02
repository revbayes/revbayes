#ifndef Func_TamuraNei_H
#define Func_TamuraNei_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Tamura-Nei rate matrix function.
     *
     * The RevLanguage wrapper of the Tamura-Nei three parameter rate matrix connects
     * the variables/parameters of the function and creates the internal TamuraNeiRateMatrixFunction object.
     * Please read the TamuraNeiRateMatrixFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-14, version 1.0
     *
     */
    class Func_TamuraNei : public TypedFunction<RateMatrix> {
        
    public:
        Func_TamuraNei( void );
        
        // Basic utility functions
        Func_TamuraNei*                                                     clone(void) const;                                          //!< Clone the object
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
