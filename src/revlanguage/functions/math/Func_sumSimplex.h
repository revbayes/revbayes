#ifndef Func_sumSimplex_H
#define Func_sumSimplex_H

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the sum function for positive real numbers.
     *
     * The RevLanguage wrapper of the sum function connects
     * the variables/parameters of the function and creates the internal SumFunction object.
     * Please read the SFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-13, version 1.0
     *
     */
    class Func_sumSimplex :  public TypedFunction<RealPos> {
        
    public:
        Func_sumSimplex( void );
        
        // Basic utility functions
        Func_sumSimplex*                                clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
