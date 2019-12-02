#ifndef Func_pruneTree_h
#define Func_pruneTree_h

#include "RlTree.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_pruneTree :  public TypedFunction<Tree> {
        
    public:
        Func_pruneTree( void );
        
        // Basic utility functions
        Func_pruneTree*                                                  clone(void) const;                                          //!< Clone the object
        virtual RevPtr<RevVariable>                                         execute(void);
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::Tree>*                    createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif /* Func_pruneTree_h */
