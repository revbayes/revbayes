#ifndef Func_powerVector_H
#define Func_powerVector_H

#include <iosfwd>
#include <vector>

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedFunction.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    
    /**
     * The RevLanguage wrapper of the exponentiation/power function for vectors.
     *
     * The RevLanguage wrapper of the power function connects
     * the variables/parameters of the function and creates the internal PowerVectorFunction object.
     * Please read the PowerVectorFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-29, version 1.0
     *
     */
    class Func_powerVector : public TypedFunction< ModelVector< RealPos > > {
        
    public:
        Func_powerVector( void );
        
        // Basic utility functions
        Func_powerVector*                                                           clone(void) const;                                          //!< Clone the object
        static const std::string&                                                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RbVector< double > >*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
