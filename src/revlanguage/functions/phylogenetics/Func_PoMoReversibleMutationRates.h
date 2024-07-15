/**
 * This file contains the declaration of the RevLanguage PoMoReversibleMutationRates function, which
 * is used to create a vector of root frequencies for the PoMo model.
 */


#ifndef Func_PoMoReversibleMutationRates_H
#define Func_PoMoReversibleMutationRates_H

#include <string>
#include <iosfwd>
#include <vector>

#include "RlSimplex.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "ModelVector.h"


namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_PoMoReversibleMutationRates : public TypedFunction< ModelVector<RealPos> > {
        
    public:
        Func_PoMoReversibleMutationRates( void );
        
        // Basic utility functions
        Func_PoMoReversibleMutationRates*                                       clone(void) const;                                          //!< Clone the object
        static const std::string&                                               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                         getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >*          createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                    getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
