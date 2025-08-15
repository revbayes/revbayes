#ifndef Func_PseudoDataSub2_H
#define Func_PseudoDataSub2_H

#include <string>
#include <iosfwd>
#include <vector>

#include "Real.h"
#include "RlPseudoData.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;

    /**
     * The RevLanguage wrapper of Sub2 pseudo-data.
     *
     * @copyright Copyright 2025-
     * @author Benjamin D. Redelings
     * @since 2025-01-17, version 1.0
     *
     */
    class Func_PseudoDataSub2: public TypedFunction<PseudoData<Real>> {

    public:
        Func_PseudoDataSub2( void );

        // Basic utility functions
        Func_PseudoDataSub2*                                                clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::PseudoData<double>>*     createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#endif


