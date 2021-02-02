/**
 * @file
 * This file contains the declaration of the RevLanguage polymorphicRootFrequencies function, which
 * is used to create a vector of root frequencies for the PoMo model.
 *
 * @brief Declaration and implementation of polymorphicRootFrequencies
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_pomoRootFrequencies.h 1406 2014-08-22 02:06:14Z Boussau $
 */


#ifndef Func_pomoRootFrequencies_H
#define Func_pomoRootFrequencies_H

#include <string>
#include <iosfwd>
#include <vector>

#include "RlSimplex.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_pomoRootFrequencies : public TypedFunction<Simplex> {
        
    public:
        Func_pomoRootFrequencies( void );
        
        // Basic utility functions
        Func_pomoRootFrequencies*                                               clone(void) const;                                          //!< Clone the object
        static const std::string&                                               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                         getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::Simplex >*                   createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                    getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
