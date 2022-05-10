#ifndef Func_pomoStateConverter_H
#define Func_pomoStateConverter_H

#include <iosfwd>

#include "Procedure.h"
#include "RevPtr.h"

namespace RevLanguage {

class ArgumentRules;
class RevVariable;
class TypeSpec;
    
    class Func_pomoStateConverter : public Procedure {
        
    public:
        Func_pomoStateConverter( void );
        
        // Basic utility functions
        Func_pomoStateConverter*                        clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevPtr<RevVariable>                             execute(void);                                              //!< Execute function
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                  //!< Get type of return value
        
    };
    
}

#endif
