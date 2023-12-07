#ifndef Func_sigmoidVector_H
#define Func_sigmoidVector_H

#include "RealPos.h"
#include "RlTypedFunction.h"
#include "ModelVector.h"

#include <string>

namespace RevLanguage {
    
    class Func_sigmoidVector :  public TypedFunction< ModelVector< RealPos > > {
        
    public:
        Func_sigmoidVector( void );
        
        // Basic utility functions
        Func_sigmoidVector*                            clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                          getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
