/**
 * @file
 * This file contains the declaration of Func_module, which is
 * the function used to read commands (module) from a file.
 *
 * @brief Declaration of Func_module
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-22 12:06:58 -0700 (Sun, 22 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @interface RbFunction
 * @package functions
 * @since Version 1.0, 2009-09-03
 *
 * $Id: Func_module.h 1431 2012-04-22 19:06:58Z hoehna $
 */

#ifndef Func_module_H
#define Func_module_H

#include "Procedure.h"

namespace RevLanguage {
    
    class Func_module : public Procedure {
        
    public:
        Func_module( void );
        
        // Basic utility functions
        Func_module*                            clone(void) const;                                          //!< Clone object
        static const std::string&               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                         getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Func_module functions
        const ArgumentRules&                    getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                         getReturnType(void) const;                                  //!< Get type of return val
        
        RevPtr<RevVariable>                     execute(void);                                              //!< Execute function
        
    };
    
}

#endif

