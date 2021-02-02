/**
 * @file
 * This file contains the declaration of the RevLanguage average distance matrix
 * function, which is used to create a deterministic variable associated with the
 * average distance matrix.
 *
 * @brief Declaration and implementation of Func_gtr
 *
 * @author David Cerny
 * @license GPL version 3
 * @version 1.0
 * @since 2019-10-23, version 1.0
 *
 */


#ifndef Func_avgDistanceMatrix_H
#define Func_avgDistanceMatrix_H

#include <string>
#include <iosfwd>
#include <vector>

#include "RlAverageDistanceMatrix.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "AverageDistanceMatrix.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_avgDistanceMatrix :  public TypedFunction<AverageDistanceMatrix> {
        
    public:
        Func_avgDistanceMatrix( void );
        
        // Basic utility functions
        Func_avgDistanceMatrix*                                             clone(void) const;                     //!< Clone the object
        static const std::string&                                           getClassType(void);                    //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                //!< Get class type spec
        std::vector<std::string>                                            getFunctionNameAliases(void) const;    //!< Get aliases for the Rev name of the function
        std::string                                                         getFunctionName(void) const;           //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;               //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::AverageDistanceMatrix>*   createFunction(void) const;            //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;          //!< Get argument rules
        
    };
    
}

#endif
