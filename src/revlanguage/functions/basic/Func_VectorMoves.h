#ifndef Func_VectorMoves_H
#define Func_VectorMoves_H

#include "Procedure.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the VectorMoves function.
     *
     * The RevLanguage wrapper of the VectorMoves function connects
     * the variables/parameters of the function and creates the internal VectorMovesFunction object.
     * Please read the VectorMovesFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_VectorMoves : public Procedure {
        
    public:
        Func_VectorMoves();
        
        // Basic utility functions
        Func_VectorMoves*                                                           clone(void) const;                                          //!< Clone the object
        virtual RevPtr<RevVariable>                                                 execute(void);                                              //!< Create a random variable from this distribution
        const ArgumentRules&                                                        getArgumentRules(void) const;                               //!< Get argument rules
        static const std::string&                                                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getReturnType(void) const;                                                      //!< Get type of return value
        const TypeSpec&                                                             getTypeSpec(void) const;                                    //!< Get language type of the object

        
    };
    
}


#endif


