/**
 * @file
 * This file contains the declaration of Func_readAncestralStateTreeTrace.
 * This class is differentiated from Func_readTreeTrace only in the sense that
 * the trees put into the TraceTree have not been rerooted; the nodes keep the
 * indexes already assigned to them, which is necessary for ancestral state
 * reconstruction.
 *
 * @brief Declaration of Func_readAncestralStateTreeTrace
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 */


#ifndef Func_readAncestralStateTreeTrace_H
#define Func_readAncestralStateTreeTrace_H

#include "RlBranchLengthTree.h"
#include "Procedure.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {
    
    class Func_readAncestralStateTreeTrace :  public Procedure {
        
    public:
        // Basic utility functions
        Func_readAncestralStateTreeTrace*       clone(void) const;                                                                  //!< Clone the object
        static const std::string&               getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                  getClassTypeSpec(void);                                                             //!< Get class type spec
        std::string                             getFunctionName(void) const;                                                        //!< Get the primary name of the function in Rev
        const TypeSpec&                         getTypeSpec(void) const;                                                            //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                     execute(void);                                                                      //!< Execute function
        const ArgumentRules&                    getArgumentRules(void) const;                                                       //!< Get argument rules
        const TypeSpec&                         getReturnType(void) const;                                                          //!< Get type of return value
        
    private:
		
        TraceTree*                              readBranchLengthTrees(const std::vector<RevBayesCore::path> &fns, const std::string &d);
        TraceTree*                              readTimeTrees(const std::vector<RevBayesCore::path> &fns, const std::string &d);
    };
    
}

#endif

