/**
 * @file
 * This file contains the declaration of the RevLanguage minimum branch length
 * time-scaling function, which time-scales a topology based on a vector of tip
 * ages and a user-specified minimum branch length.
 *
 * @brief Declaration and implementation of Func_MinBLTimeScaling
 *
 * @author David Cerny, Laura Mulvey
 * @license GPL version 3
 * @version 1.2.6
 * @since 2025-01-09, version 1.2.6
 *
 */

#ifndef Func_MinBLTimeScaling_H
#define Func_MinBLTimeScaling_H

#include <string>

#include "Procedure.h"
#include "RevPtr.h"
#include "Tree.h"
#include "TopologyNode.h"

namespace RevLanguage {
    
    class Func_MinBLTimeScaling :  public Procedure {
        
    public:
        Func_MinBLTimeScaling( void );
        
        // Basic utility functions
        Func_MinBLTimeScaling*          clone(void) const;                     //!< Clone the object
        static const std::string&       getClassType(void);                    //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                //!< Get class type spec
        std::string                     getFunctionName(void) const;           //!< Get the primary name of the function in Rev
        const TypeSpec&                 getTypeSpec(void) const;               //!< Get the type spec of the instance
        
        // Regular functions
        RevPtr<RevVariable>             execute(void);                         //!< Execute function
        const ArgumentRules&            getArgumentRules(void) const;          //!< Get argument rules
        const TypeSpec&                 getReturnType(void) const;             //!< Get type of return value
    };
    
}

#endif
