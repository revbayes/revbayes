#ifndef Func_writePhylip_H
#define Func_writePhylip_H

#include "Procedure.h"


namespace RevLanguage {

    /**
     * Function that takes in character data and writes it into a file in Phylip format.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2023-06-05, version 1.2.1
     */
    class Func_writePhylip : public Procedure {
        
    public:
        // Basic utility functions
        Func_writePhylip*                   clone(void) const;                                          //!< Clone the object
        static const std::string&           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                 execute(void);                                              //!< Execute function
        const ArgumentRules&                getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                     getReturnType(void) const;                                  //!< Get type of return value
        
        
    };
    
}

#endif

