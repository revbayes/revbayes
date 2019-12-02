#ifndef Func_simplex_H
#define Func_simplex_H

#include "RlSimplex.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the simplex function.
     *
     * The RevLanguage wrapper of the simplex function connects
     * the variables/parameters of the function and creates the internal SimplexFunction object.
     * Please read the SimplexFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_simplex : public TypedFunction<Simplex> {
        
    public:
        Func_simplex();
        
        // Basic utility functions
        Func_simplex*                                                               clone(void) const;                                          //!< Clone the object
        static const std::string&                                                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::Simplex >*                       createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif


