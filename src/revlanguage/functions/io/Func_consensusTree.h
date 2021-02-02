/**
 * @file
 * This file contains the declaration of Func_consensusTree, which calculates the majority-rule consensus tree from a tree
 * trace and adds the node posterior probabilities.
 *
 * @brief Declaration of Func_consensusTree
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-27 16:10:36 +0200 (Sun, 27 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_consensusTree.h 1587 2012-05-27 14:10:36Z seraklop $
 */

#ifndef Func_consensusTree_H
#define Func_consensusTree_H

#include "Procedure.h"
#include "RbFileManager.h"

#include <string>


namespace RevLanguage {

    class Func_consensusTree :  public Procedure {
        
    public:
        // Basic utility functions
		Func_consensusTree*                 clone(void) const;                                          //!< Clone the object
        static const std::string&           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        std::vector<std::string>            getFunctionNameAliases( void ) const;                       //!< Get aliases for the Rev name of the function
        const TypeSpec&                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                 execute(void);                                              //!< Execute function
        const ArgumentRules&                getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                     getReturnType(void) const;                                  //!< Get type of return value
        
        
    };
    
}

#endif

