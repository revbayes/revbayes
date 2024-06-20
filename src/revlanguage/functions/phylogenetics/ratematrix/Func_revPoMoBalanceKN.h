#ifndef Func_revPoMoBalanceKN_H
#define Func_revPoMoBalanceKN_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the revPoMoBalanceKN rate matrix function.
     *
     * The RevLanguage wrapper of the revPoMoBalanceKN rate matrix connects
     * the variables/parameters of the function and creates the internal Func_revPoMoBalanceKN object.
     *
     * @brief Declaration of RateMatrix_revPoMoBalanceKN, a reversible matrix combining polymorphisms,
     * substitutions and the balancing selection
     *
     * @copyright Copyright 2023-
     * @author The RevBayes Development Core Team (Svitlana Braichenko)
     * @since 2023-12-12, version 1.2.2
     *
     */
    class Func_revPoMoBalanceKN : public TypedFunction<RateMatrix> {
        
    public:
        Func_revPoMoBalanceKN( void );
        
        // Basic utility functions
        Func_revPoMoBalanceKN*                                                   clone(void) const;                                          //!< Clone the object
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


