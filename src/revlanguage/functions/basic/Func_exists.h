#ifndef Func_exists_H
#define Func_exists_H

#include "Procedure.h"

namespace RevLanguage {
    
    /**
     * @brief Rev function to check if a variable exists in the workspace.
     *
     * Check in the workspace if a variable exists and returns either true or false.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-10-09
     *
     */
    class Func_exists : public Procedure {
        
    public:
        Func_exists( void );
        
        // Basic utility functions
        Func_exists*                                    clone(void) const;                                                          //!< Clone object
        static const std::string&                       getClassType(void);                                                         //!< Get Rev exists
        static const TypeSpec&                          getClassTypeSpec(void);                                                     //!< Get class exists spec
        std::string                                     getFunctionName(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                    //!< Get language exists of the object
        
        // Func_source functions
        const ArgumentRules&                            getArgumentRules(void) const;                                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                                  //!< Get exists of return val
        
        RevPtr<RevVariable>                             execute(void);                                                              //!< Execute function
        
    protected:
        
    };
    
}

#endif

