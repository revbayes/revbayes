#ifndef Func_readTreeTrace_H
#define Func_readTreeTrace_H

#include "Clade.h"
#include "Procedure.h"
#include "RlBranchLengthTree.h"
#include "RlClade.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"
#include "RbFileManager.h"
#include "WorkspaceVector.h"

#include <map>
#include <cstdint>
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
        
        WorkspaceVector<TraceTree>*         readTrees(const std::vector<RevBayesCore::path> &fns, const std::string &d, const std::string& treetype, bool unroot_nonclock, std::int64_t thin, std::int64_t offset);
        WorkspaceVector<TraceTree>*         readTreesNexus(const std::vector<RevBayesCore::path> &fns, const std::string& treetype, bool unroot_nonclock, std::int64_t thin, std::int64_t offset);  //!< Read tree trace from Nexus file(s)
    };
    
}

#endif

