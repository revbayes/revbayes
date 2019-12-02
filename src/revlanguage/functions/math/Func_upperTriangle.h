#ifndef Func_upperTriangle_H
#define Func_upperTriangle_H

#include "RlMatrixRealSymmetric.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_upperTriangle : public TypedFunction< ModelVector< RealPos > > {
        
    public:
        Func_upperTriangle();
        
        // Basic utility functions
        Func_upperTriangle*                                                 clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassName(void);                                         //!< Get class name
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >*      createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif


