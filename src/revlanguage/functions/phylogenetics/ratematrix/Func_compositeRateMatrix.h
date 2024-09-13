//
//  Func_compositeRateMatrix.h
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/13/24.
//

#ifndef Func_compositeRateMatrix_h
#define Func_compositeRateMatrix_h

#include <string>
#include <iosfwd>
#include <vector>

#include "RateGenerator.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "RlRateGenerator.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_compositeRateMatrix : public TypedFunction<RateGenerator> {
        
    public:
        Func_compositeRateMatrix( void );
        
        // Basic utility functions
        Func_compositeRateMatrix*                                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >*     createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif /* Func_compositeRateMatrix_h */
