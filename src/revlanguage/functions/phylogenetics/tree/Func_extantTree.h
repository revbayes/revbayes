#ifndef Func_extantTree_h
#define Func_extantTree_h

#include "RlTimeTree.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_extantTree :  public TypedFunction<TimeTree> {
        
    public:
        Func_extantTree( void );
        
        // Basic utility functions
        Func_extantTree*                                                    clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::Tree>*                    createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif /* Func_extantTree_h */
