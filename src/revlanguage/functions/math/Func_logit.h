#ifndef Func_logit_H
#define Func_logit_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the logit function.
     *
     * @copyright Copyright 2024-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2024-07-11, version 1.0
     *
     */
    class Func_logit : public TypedFunction<Real> {

    public:
        Func_logit( void );

        // Basic utility functions
        Func_logit*                                     clone(void) const;                                          //!< Clone the object
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
