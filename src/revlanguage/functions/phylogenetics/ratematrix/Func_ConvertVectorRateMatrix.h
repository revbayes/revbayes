#ifndef Func_ConvertVectorRateMatrix_H
#define Func_ConvertVectorRateMatrix_H

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

#include "ModelVector.h"
#include "RbVector.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;

    /**
     * RevLanguage type conversion function from RateMatrix[] -> SiteMixtureModel[]
     *
     * @copyright Copyright 2024-
     * @author Benjamin D. Redelings
     * @since 2024-04-17, version 1.0
     *
     */
    class Func_ConvertVectorRateMatrix: public TypedFunction< ModelVector<SiteMixtureModel> > {

    public:
        Func_ConvertVectorRateMatrix( void );

        // Basic utility functions
        Func_ConvertVectorRateMatrix*                                       clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<RevBayesCore::SiteMixtureModel> >*          createFunction(void) const;                                 //!< Create a function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#endif


