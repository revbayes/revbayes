//
//  Func_TreeTrace.h
//  RevBayes_development_branch
//
//  Created by will freyman on 12/15/15.
//  Copyright (c) 2015 will freyman. All rights reserved.
//

#ifndef Func_treeTrace_H
#define Func_treeTrace_H

#include "RlBranchLengthTree.h"
#include "Procedure.h"
#include "RbFileManager.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {
    
    class Func_treeTrace : public Procedure {
        
    public:
        // Basic utility functions
        Func_treeTrace*                     clone(void) const;                                                                  //!< Clone the object
        static const std::string&           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                             //!< Get class type spec
        std::string                         getFunctionName(void) const;                                                        //!< Get the primary name of the function in Rev
        const TypeSpec&                     getTypeSpec(void) const;                                                            //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                 execute(void);                                                                      //!< Execute function
        const ArgumentRules&                getArgumentRules(void) const;                                                       //!< Get argument rules
        std::vector<std::string>            getHelpAuthor(void) const;      
        std::string                         getHelpDescription(void) const; 
        std::string                         getHelpExample(void) const;    
        std::vector<std::string>            getHelpSeeAlso(void) const;     
        const TypeSpec&                     getReturnType(void) const;                                                          //!< Get type of return value
        
    private:

    };
    
}

#endif
