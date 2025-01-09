/**
 * @file
 * This file contains the declaration of the RevLanguage minimum branch length
 * time-scaling function, which is used to time-scale a topology based on the ages
 * of its terminal taxa.
 *
 * @brief Declaration and implementation of Func_gtr
 *
 * @author David Cerny, Laura Mulvey
 * @license GPL version 3
 * @version 1.0
 * @since 2025-01-09, version 1.0
 *
 */


#ifndef Func_MinBLTimeScaling_H
#define Func_MinBLTimeScaling_H

#include <string>
#include <iosfwd>
#include <vector>

#include "RlTree.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "MinBLTimeScalingFunction.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "Tree.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_MinBLTimeScaling :  public TypedFunction<Tree> {
        
    public:
        Func_MinBLTimeScaling( void );
        
        // Basic utility functions
        Func_MinBLTimeScaling*                                              clone(void) const;                     //!< Clone the object
        static const std::string&                                           getClassType(void);                    //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                //!< Get class type spec
        std::string                                                         getFunctionName(void) const;           //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;               //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::Tree>*                    createFunction(void) const;            //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;          //!< Get argument rules
        
    };
    
}

#endif
