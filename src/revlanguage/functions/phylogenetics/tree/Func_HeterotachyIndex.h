/**
 * @file
 * This file contains the declaration of the RevLanguage heterotachyIndex statistic, which
 * is used to create a deterministic variable associated with the heterotachyIndex statistic.
 *
 * @brief Declaration of the Heterotachy Index of Wang, Susko and Roger (2019)
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2023-10-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_heterotachyIndex.h 1406 2023-10-20 02:06:14Z Wright, Mulvey, Khakurel $
 */


#ifndef Func_heterotachyIndex_H
#define Func_heterotachyIndex_H

#include <string>
#include <iosfwd>
#include <vector>

#include "RealPos.h"
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
    
    class Func_heterotachyIndex : public TypedFunction<RealPos> {
        
    public:
        Func_heterotachyIndex( void );
        
        // Basic utility functions
        Func_heterotachyIndex*                       clone(void) const;                                          //!< Clone the object
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
