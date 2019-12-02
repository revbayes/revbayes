#ifndef Func_listOptions_H
#define Func_listOptions_H

#include "Procedure.h"

namespace RevLanguage {
    
    
    /**
     * @brief Rev function to get a user option of RevBayes.
     *
     * The listOptions function lets the user get/see a list of all options.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-02-01
     *
     */
    class Func_listOptions : public Procedure {
        
    public:
        Func_listOptions( void );
        
        // Basic utility functions
        Func_listOptions*                               clone(void) const;                                                          //!< Clone object
        static const std::string&                       getClassType(void);                                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                    //!< Get language type of the object
        
        // Func_source functions
        const ArgumentRules&                            getArgumentRules(void) const;                                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                                  //!< Get type of return val
        
        RevPtr<RevVariable>                             execute(void);                                                              //!< Execute function
        
    protected:
        
    };
    
}

#endif

