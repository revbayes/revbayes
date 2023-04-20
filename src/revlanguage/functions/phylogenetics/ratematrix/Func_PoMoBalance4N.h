#ifndef Func_PoMoBalance4N_H
#define Func_PoMoBalance4N_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the PoMoBalance4N rate matrix function.
     *
     * The RevLanguage wrapper of the PoMoBalance4N rate matrix connects
     * the variables/parameters of the function and creates the internal Func_PoMoBalance4N object.
     *
     * @brief Declaration of RateMatrix_PoMoBalance4N, a reversible matrix combining polymorphisms, substitutions and the balancing selection
     *
     * @copyright Copyright 2022-
     * @author The RevBayes Development Core Team
     * @since 2022-08-01, version 1.0
     *
     */
    class Func_PoMoBalance4N : public TypedFunction<RateMatrix> {
        
    public:
        Func_PoMoBalance4N( void );
        
        // Basic utility functions
        Func_PoMoBalance4N*                                                        clone(void) const;                                          //!< Clone the object
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


