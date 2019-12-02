#ifndef Func_simplexFromVector_H
#define Func_simplexFromVector_H

#include "RlTypedFunction.h"
#include "RlSimplex.h"

#include <string>

namespace RevLanguage {
    
    class Func_simplexFromVector :  public TypedFunction<Simplex> {
        
    public:
        Func_simplexFromVector();
        
        // Basic utility functions
        Func_simplexFromVector*                                                 clone(void) const;                                          //!< Clone the object
        static const std::string&                                               getClassType(void);                                         //!< Get Rev type (static)
        static const TypeSpec&                                                  getClassTypeSpec(void);                                     //!< Get Rev type spec (static)
        std::string                                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                         getTypeSpec(void) const;                                    //!< Get Rev type spec (from instance)
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::Simplex >*                   createFunction(void) const;                                                     //!< Create internal function object
        const ArgumentRules&                                                    getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif


