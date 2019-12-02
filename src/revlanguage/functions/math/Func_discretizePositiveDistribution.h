#ifndef Func_discretizePositiveDistribution_H
#define Func_discretizePositiveDistribution_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_discretizePositiveDistribution :  public TypedFunction< ModelVector< RealPos > > {
        
    public:
        Func_discretizePositiveDistribution();
        
        // Basic utility functions
        Func_discretizePositiveDistribution*                                clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassName(void);                                         //!< Get class name
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >*      createFunction(void) const;                                                     //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif


