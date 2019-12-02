#ifndef Func_mapTree_H
#define Func_mapTree_H

#include "Procedure.h"
#include "RbFileManager.h"

#include <string>


namespace RevLanguage {
    
    class Func_mapTree : public Procedure {
        
    public:
        // Basic utility functions
        Func_mapTree*                       clone(void) const;                                          //!< Clone the object
        static const std::string&           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                     //!< Get class type spec
        std::vector<std::string>            getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                         getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                         getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<std::string>            getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
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

