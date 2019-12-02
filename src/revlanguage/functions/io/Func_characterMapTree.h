//
//  Func_characterMapTree.h
//  RevBayes_development_branch
//
//  Created by will freyman on 12/22/16.
//  Copyright (c) 2016 will freyman. All rights reserved.
//

#ifndef __RevBayes_development_branch__Func_characterMapTree__
#define __RevBayes_development_branch__Func_characterMapTree__

#include "Procedure.h"
#include "RbFileManager.h"

#include <string>


namespace RevLanguage {
    
    class Func_characterMapTree : public Procedure {
        
    public:
        // Basic utility functions
        Func_characterMapTree*              clone(void) const;                                          //!< Clone the object
        static const std::string&           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                 execute(void);                                              //!< Execute function
        const ArgumentRules&                getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                     getReturnType(void) const;                                  //!< Get type of return value
        
    private:
        
    };
    
}

#endif
