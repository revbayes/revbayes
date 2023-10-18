#ifndef Func_readAtlas_H
#define Func_readAtlas_H

#include "Procedure.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {
    
    class Func_readAtlas :  public Procedure {
        
    public:
        // Basic utility functions
        Func_readAtlas*                     clone(void) const;                                          //!< Clone the object
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

