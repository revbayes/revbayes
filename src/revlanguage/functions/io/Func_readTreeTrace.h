/**
 * @file
 * This file contains the declaration of Func_readTreeTrace, which uses the NCL to read in trees in any of the supported formats:
 * nexus
 *
 * @brief Declaration of Func_readTrees
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-27 16:10:36 +0200 (Sun, 27 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_readDiscreteCharacterData.h 1587 2012-05-27 14:10:36Z hoehna $
 */

#ifndef Func_readTreeTrace_H
#define Func_readTreeTrace_H

#include "Clade.h"
#include "Procedure.h"
#include "RlBranchLengthTree.h"
#include "RlClade.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {
    
    class Func_readTreeTrace : public Procedure {
        
    public:
        // Basic utility functions
        Func_readTreeTrace*                 clone(void) const;                                                                  //!< Clone the object
        static const std::string&           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                             //!< Get class type spec
        std::string                         getFunctionName(void) const;                                                        //!< Get the primary name of the function in Rev
        const TypeSpec&                     getTypeSpec(void) const;                                                            //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                 execute(void);                                                                      //!< Execute function
        const ArgumentRules&                getArgumentRules(void) const;                                                       //!< Get argument rules
        const TypeSpec&                     getReturnType(void) const;                                                          //!< Get type of return value
        
    private:
        
        TraceTree*                          readTrees(const std::vector<std::string> &fns, const std::string &d, bool clock, long thin);
        TraceTree*                          readTreesNexus(const std::vector<std::string> &fns, bool clock, long thin);  //!< Read tree trace from Nexus file(s)
    };
    
}

#endif

