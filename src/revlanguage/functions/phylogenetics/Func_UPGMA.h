#ifndef Func_UPGMA_H
#define Func_UPGMA_H

#include "Procedure.h"

namespace RevLanguage {
    
    /**
     * @brief Rev function to construct a UPGMA tree..
     *
     * This procedure builds a tree using the UPGMA algorithm.
     *
     *
     * @copyright Copyright 2009-
     * @author Sebastian Hoehna
     * @since Version 1.0, 2014-12-09
     *
     */
    class Func_UPGMA : public Procedure {
        
    public:
        Func_UPGMA( void );
        
        // Basic utility functions
        Func_UPGMA*                                     clone(void) const;                                          //!< Clone object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Func_source functions
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                  //!< Get type of return val
        
        RevPtr<RevVariable>                             execute(void);                                              //!< Execute function

    protected:

    };
    
}

#endif

