#ifndef Func_PoMoKN_H
#define Func_PoMoKN_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the PoMo-KN rate matrix function.
     *
     * The RevLanguage wrapper of the PoMo-KN rate matrix connects
     * the variables/parameters of the function and creates the internal PoMoKNRateMatrixFunction object.
     * Please read the PoMoKNRateMatrixFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-14, version 1.0
     *
     */
    class Func_PoMoKN : public TypedFunction<RateMatrix> {
        
    public:
        Func_PoMoKN( void );
        
        // Basic utility functions
        Func_PoMoKN*                                                        clone(void) const;                                          //!< Clone the object
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


