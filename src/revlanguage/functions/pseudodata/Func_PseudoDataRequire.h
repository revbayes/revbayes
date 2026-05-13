#ifndef Func_PseudoDataRequire_H
#define Func_PseudoDataRequire_H

#include <string>
#include <iosfwd>
#include <vector>

#include "Real.h"
#include "RlPseudoDataLikelihood.h"
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
     * The RevLanguage wrapper of Require pseudo-data.
     *
     * @copyright Copyright 2026-
     * @author Benjamin D. Redelings
     * @since 2026-0-27, version 1.0
     *
     */
    class Func_PseudoDataLikelihoodRequire: public TypedFunction<PseudoDataLikelihood> {

    public:
        Func_PseudoDataLikelihoodRequire( void ): TypedFunction<PseudoDataLikelihood>( ) { }

        // Basic utility functions
        Func_PseudoDataLikelihoodRequire* clone( void ) const
        {
            return new Func_PseudoDataLikelihoodRequire( *this );
        }

        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::PseudoDataLikelihood>*   createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}
#endif
