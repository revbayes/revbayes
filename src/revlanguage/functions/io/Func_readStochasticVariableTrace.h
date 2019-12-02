/**
 * @file
 * This file contains the declaration of Func_readStochasticVariableTrace
 *
 * @brief Declaration of Func_readStochasticVariableTrace
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-27 16:10:36 +0200 (Sun, 27 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_readDiscreteCharacterData.h 1587 2012-05-27 14:10:36Z hoehna $
 */

#ifndef Func_readStochasticVariableTrace_H
#define Func_readStochasticVariableTrace_H

#include "Procedure.h"
#include "RbFileManager.h"
#include "RlTraceTree.h"

#include <string>
#include <vector>


namespace RevLanguage {
    
    class Func_readStochasticVariableTrace : public Procedure {
        
    public:
        // Basic utility functions
        Func_readStochasticVariableTrace*       clone(void) const;                                                      //!< Clone the object
        static const std::string&               getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&                  getClassTypeSpec(void);                                                 //!< Get class type spec
        std::string                             getFunctionName(void) const;                                            //!< Get the primary name of the function in Rev
        const TypeSpec&                         getTypeSpec(void) const;                                                //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                     execute(void);                                                          //!< Execute function
        const ArgumentRules&                    getArgumentRules(void) const;                                           //!< Get argument rules
        const TypeSpec&                         getReturnType(void) const;                                              //!< Get type of return value
        
    };
    
}

#endif

