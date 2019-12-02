#ifndef Func_normalize_H
#define Func_normalize_H

#include "RlTypedFunction.h"
#include "ModelVector.h"
#include "RealPos.h"
#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_normalize : public TypedFunction< ModelVector<RealPos> > {
        
    public:
        Func_normalize();
        
        // Basic utility functions
        Func_normalize*                                                     clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >*      createFunction(void) const;                                 //!< Create internal function object
        
    };
    
}


#endif


