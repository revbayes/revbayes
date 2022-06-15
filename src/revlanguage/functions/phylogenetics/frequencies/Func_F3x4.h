#ifndef Func_F3x4_H
#define Func_F3x4_H

#include <string>
#include <iosfwd>
#include <vector>

#include "RlSimplex.h"
#include "RlRateMatrix.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RateGenerator.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;

    /**
     * The RevLanguage wrapper for the F3x4 rate matrix function.
     *
     * The RevLanguage wrapper
     * + connects the variables/parameters of the function to the F3x4Function object
     *   (see F3x4Func.{h,cpp})
     * + the F3x4Func object then creates a Simplex object.
     *
     * @copyright Copyright 2021-
     * @author Benjamin D. Redelings
     * @since 2021-11-24, version 1.0
     *
     */
    class Func_F3x4 : public TypedFunction< Simplex > {

    public:
        // Basic utility functions
        Func_F3x4*                                                          clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::Simplex >*               createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#endif


