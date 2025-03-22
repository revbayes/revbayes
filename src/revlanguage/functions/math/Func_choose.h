//
//  Func_choose.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/21/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Func_choose_h
#define Func_choose_h

#include <string>
#include <iosfwd>
#include <vector>

#include "Natural.h"
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
     * The RevLanguage wrapper of the choose function.
     *
     * The RevLanguage wrapper of the choose function connects
     * the variables/parameters of the function and creates the internal ChooseFunction object.
     * Please read the ChooseFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2016-10-21, version 1.0
     *
     */
    class Func_choose : public TypedFunction<Natural> {
        
    public:
        Func_choose( void );
        
        // Basic utility functions
        Func_choose*                                     clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<std::int64_t>*              createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules

    };
    
}


#endif /* Func_choose_h */
