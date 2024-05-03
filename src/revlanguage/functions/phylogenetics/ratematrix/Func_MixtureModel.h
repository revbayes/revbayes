#ifndef Func_MixtureModel_H
#define Func_MixtureModel_H

#include <string>
#include <iosfwd>
#include <vector>

#include "RlSiteMixtureModel.h"
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
     * The RevLanguage wrapper of the Gamma model of rates-across-sites.
     *
     * @copyright Copyright 2021-
     * @author Benjamin D. Redelings
     * @since 2022-09-22, version 1.0
     *
     */
    class Func_MixtureModel: public TypedFunction<SiteMixtureModel> {

    public:
        Func_MixtureModel( void );

        // Basic utility functions
        Func_MixtureModel*                                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::SiteMixtureModel >*      createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#endif


