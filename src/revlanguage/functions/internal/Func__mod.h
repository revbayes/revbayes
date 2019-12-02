#ifndef Func__mod_H
#define Func__mod_H

#include "Natural.h"
#include "RlTypedFunction.h"
#include <string>

namespace RevLanguage {
    
    
    /**
     * @brief Modulo function.
     *
     * Performs the modulo operation of two natural numbers, e.g.,
     * 7 % 5 = 2
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-07-04
     *
     */
    class Func__mod: public TypedFunction<Natural> {
        
    public:
        Func__mod();
        
        // Basic utility functions
        Func__mod*                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                             getTypeSpec(void) const;                                    //!< Get language type of the object
        bool                                        isInternal(void) const { return true; }                     //!< Is this an internal function?

        // Regular functions
        RevBayesCore::TypedFunction<long>*           createFunction(void) const;                                                     //!< Create internal function object
        const ArgumentRules&                        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif


