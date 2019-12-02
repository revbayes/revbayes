//
//  Func_convertToPhylowood.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__Func_convertToPhylowood__
#define __revbayes_proj__Func_convertToPhylowood__

#include "Procedure.h"
#include "RbFileManager.h"

#include <string>


namespace RevLanguage {
    
    class Func_convertToPhylowood :  public Procedure {
        
    public:
        // Basic utility functions
        Func_convertToPhylowood*            clone(void) const;                                          //!< Clone the object
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


#endif /* defined(__revbayes_proj__Func_convertToPhylowood__) */
