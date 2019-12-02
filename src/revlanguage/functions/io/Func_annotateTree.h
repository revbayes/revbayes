/**
 * @file
 * This file contains the declaration of Func_annotateTree, which annotates the
 * x% Highest Posterior Density (HPD) interval of node ages onto the input TimeTree.
 *
 * @brief Declaration of Func_annotateTree
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 */

#ifndef Func_annotateTree_H
#define Func_annotateTree_H

#include "Procedure.h"
#include "RbFileManager.h"

#include <string>


namespace RevLanguage {
    
    class Func_annotateTree : public Procedure {
        
    public:
        // Basic utility functions
        Func_annotateTree*              clone(void) const;                                          //!< Clone the object
        static const std::string&       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                 getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>             execute(void);                                              //!< Execute function
        const ArgumentRules&            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                 getReturnType(void) const;                                  //!< Get type of return value
        
    private:
        
    };
    
}


#endif

