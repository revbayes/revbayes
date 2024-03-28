#ifndef Func_SmoothenTimeline_H
#define Func_SmoothenTimeline_H

#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RbVector.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the SmoothenTimeline function.
     *
     * The RevLanguage wrapper of the SmoothenTimeline function connects
     * the variables/parameters of the function and creates the internal vector.
     * Please read the SmoothenTimelineFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team
     * @since 2014-08-14, version 1.0
     *
     */
    class Func_SmoothenTimeline : public TypedFunction<ModelVector<RealPos> > {

    public:
        Func_SmoothenTimeline( void );

        // Basic utility functions
        Func_SmoothenTimeline*                                              clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >*      createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#endif
