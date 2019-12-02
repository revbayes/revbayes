#ifndef Func_discretizeGammaQuadrature_H
#define Func_discretizeGammaQuadrature_H

#include "RlMatrixRealPos.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_discretizeGammaQuadrature : public TypedFunction<MatrixRealPos> {
        
    public:
        Func_discretizeGammaQuadrature();
        
        // Basic utility functions
        Func_discretizeGammaQuadrature*                                     clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassName(void);                                         //!< Get class name
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif


