#ifndef Func_sigmoid_H
#define Func_sigmoid_H

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    class Func_sigmoid :  public TypedFunction<RealPos> {
        
    public:
        Func_sigmoid( void );
        
        // Basic utility functions
        Func_sigmoid*                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
