//
//  Func_matrix.h
//  revbayes-proj
//
//  Created by Michael Landis on 7/5/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef Func_matrix_h
#define Func_matrix_h

#include <string>
#include <iosfwd>
#include <vector>

#include "RlMatrixReal.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    /**
     * The RevLanguage wrapper of the simplex function.
     *
     * The RevLanguage wrapper of the simplex function connects
     * the variables/parameters of the function and creates the internal SimplexFunction object.
     * Please read the SimplexFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_matrix : public TypedFunction<MatrixReal> {
        
    public:
        Func_matrix();
        
        // Basic utility functions
        Func_matrix*                                                                clone(void) const;                                          //!< Clone the object
        static const std::string&                                                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >*                    createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif /* Func_matrix_h */
