#ifndef Func_stochasticMatrix_h
#define Func_stochasticMatrix_h

#include <string>
#include <iosfwd>
#include <vector>

#include "RlMatrixReal.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "RlStochasticMatrix.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    /**
     * The RevLanguage constructor for (row) stochastic matrices. A (row)
     * stochastic matrix has rows that sum to 1.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael R. May)
     * @since 2020-03-10, version 1.0
     *
     */
    class Func_stochasticMatrix : public TypedFunction<StochasticMatrix> {
        
    public:
        Func_stochasticMatrix();
        
        // Basic utility functions
        Func_stochasticMatrix*                                                      clone(void) const;                                          //!< Clone the object
        static const std::string&                                                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >*                    createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif /* Func_stochasticMatrix_h */
