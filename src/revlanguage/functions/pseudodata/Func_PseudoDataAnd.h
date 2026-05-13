#ifndef Func_PseudoDataAnd_H
#define Func_PseudoDataAnd_H

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
     * The RevLanguage wrapper to && pseudo-data.
     *
     * @copyright Copyright 2026-
     * @author Benjamin D. Redelings
     * @since 2026-01-37, version 1.0
     *
     */
    class Func_PseudoDataLikelihoodAnd: public TypedFunction<PseudoDataLikelihood> {

    public:
        Func_PseudoDataLikelihoodAnd( void ): TypedFunction<PseudoDataLikelihood>( ) { }

        // Basic utility functions
        Func_PseudoDataLikelihoodAnd* clone( void ) const
        {
            return new Func_PseudoDataLikelihoodAnd( *this );
        }

        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::PseudoDataLikelihood >*  createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#endif
