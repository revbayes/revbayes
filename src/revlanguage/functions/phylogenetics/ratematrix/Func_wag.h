#ifndef Func_wag_H
#define Func_wag_H

#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the WAG rate matrix function.
     *
     * The RevLanguage wrapper of the WAG rate matrix connects
     * the variables/parameters of the function and creates the internal WagRateMatrixFunction object.
     * Please read the WagRateMatrixFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-14, version 1.0
     *
     */
    class Func_wag : public TypedFunction<RateMatrix> {
        
    public:
        Func_wag( void );
        
        // Basic utility functions
        Func_wag*                                                       clone(void) const;                                          //!< Clone the object
        static const std::string&                                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RateGenerator>*      createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                               //!< Get argument rules
    private:
        Simplex* bf_empirical;
        
    };
    
}

#endif
