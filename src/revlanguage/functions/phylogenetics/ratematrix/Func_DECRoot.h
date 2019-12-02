#ifndef Func_DECRoot_H
#define Func_DECRoot_H


#include "RlTypedFunction.h"
#include "RbVector.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlSimplex.h"

#include <map>
#include <string>
#include <vector>

namespace RevLanguage {
    
    class Func_DECRoot :  public TypedFunction<Simplex> {
        
    public:
        Func_DECRoot( void );
        
        // Basic utility functions
        Func_DECRoot*                                                   clone(void) const;                                          //!< Clone the object
        static const std::string&                                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::Simplex >*           createFunction(void) const;                                 //!< Execute function
        const ArgumentRules&                                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
