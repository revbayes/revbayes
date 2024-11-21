#ifndef Func_sumNatural_H
#define Func_sumNatural_H

#include "Natural.h"
#include "RlTypedFunction.h"

#include <string>
#include <cstdint>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the sum function.
     *
     * The RevLanguage wrapper of the sum function connects
     * the variables/parameters of the function and creates the internal sumFunction object.
     * Please read the sumFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_sumNatural :  public TypedFunction<Natural> {
        
    public:
        Func_sumNatural( void );
        
        // Basic utility functions
        Func_sumNatural*                                clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<std::int64_t>*              createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
