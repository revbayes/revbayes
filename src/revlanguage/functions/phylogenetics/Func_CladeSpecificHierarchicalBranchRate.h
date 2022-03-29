#ifndef Func_CladeSpecificHierarchicalBranchRate_H
#define Func_CladeSpecificHierarchicalBranchRate_H

#include "ModelVector.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_CladeSpecificHierarchicalBranchRate :  public TypedFunction< ModelVector<RealPos> > {
        
    public:
        Func_CladeSpecificHierarchicalBranchRate( void );
        
        // Basic utility functions
        Func_CladeSpecificHierarchicalBranchRate*                           clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >*       createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
