#ifndef Func_NeighborJoining_H
#define Func_NeighborJoining_H

#include "Procedure.h"

namespace RevLanguage {
    
    /**
     * @brief Rev function to construct a NeighborJoining tree.
     *
     * This procedure builds a tree using the NeighborJoining algorithm.
     *
     *
     * @copyright Copyright 2009-
     * @author Sebastian Hoehna
     * @since Version 1.0, 2014-12-09
     *
     */
    class Func_NeighborJoining : public Procedure {
        
    public:
        Func_NeighborJoining( void );
        
        // Basic utility functions
        Func_NeighborJoining*                           clone(void) const;                                          //!< Clone object
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

