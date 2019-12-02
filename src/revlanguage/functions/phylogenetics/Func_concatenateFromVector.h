#ifndef Func_concatenateFromVector_H
#define Func_concatenateFromVector_H

#include "Procedure.h"

namespace RevLanguage {
    
    /**
     * @brief Rev function to concatenate two or more data matrices.
     *
     * This procedure concatenates the sequences of two or more data matrices.
     * The concatenation builds a new data matrix by concatenating sequences.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-12-09
     *
     */
    class Func_concatenateFromVector : public Procedure {
        
    public:
        Func_concatenateFromVector( void );
        
        // Basic utility functions
        Func_concatenateFromVector*                 clone(void) const;                                          //!< Clone object
        static const std::string&                   getClassType(void);                                         //!< Get Rev concatenate
        static const TypeSpec&                      getClassTypeSpec(void);                                     //!< Get class concatenate spec
        std::string                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                             getTypeSpec(void) const;                                    //!< Get language concatenate of the object
        
        // Func_source functions
        const ArgumentRules&                        getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                             getReturnType(void) const;                                  //!< Get concatenate of return val
        
        RevPtr<RevVariable>                         execute(void);                                              //!< Execute function
        
    };
    
}

#endif

