#ifndef Func_PhylogeneticIndependentContrasts_H
#define Func_PhylogeneticIndependentContrasts_H

#include "ModelVector.h"
#include "RbVector.h"
#include "Real.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_PhylogeneticIndependentContrasts :  public TypedFunction<ModelVector<Real> > {
        
    public:
        Func_PhylogeneticIndependentContrasts( void );
        
        // Basic utility functions
        Func_PhylogeneticIndependentContrasts*                              clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        std::vector<std::string>                                            getFunctionNameAliases( void ) const;
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >*       createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
