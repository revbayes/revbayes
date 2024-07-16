#ifndef Func_PoMoKNrecurrentMutations_H
#define Func_PoMoKNrecurrentMutations_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the PoMo-KN rate matrix function.
     *
     * The RevLanguage wrapper of the PoMo-KN rate matrix connects
     * the variables/parameters of the function and creates the internal PoMoKNrecurrentMutationsRateMatrixFunction object.
     * Please read the PoMoKNrecurrentMutationsRateMatrixFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author Rui Borges
     * @since 2024
     *
     */
    class Func_PoMoKNrecurrentMutations : public TypedFunction<RateMatrix> {
        
    public:
        Func_PoMoKNrecurrentMutations( void );
        
        // Basic utility functions
        Func_PoMoKNrecurrentMutations*                                      clone(void) const;                                          //!< Clone the object
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


