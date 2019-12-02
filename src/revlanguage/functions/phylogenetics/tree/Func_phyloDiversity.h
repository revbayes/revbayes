/**
 * @file
 * This file contains the declaration and implementation of the RevLanguage phyloDiversity function
 * which calculates the sum of all branches on the tree needed to connect the sampled taxa.
 *
 * @brief Declaration and implementation of Func_phyloDiversity
 *
 * (c) Copyright 2009- under GPL version 3
 * @author Will Freyman
 * @license GPL version 3
 * @version 1.0
 *
 */


#ifndef Func_phyloDiversity_H
#define Func_phyloDiversity_H

#include "RlTree.h"
#include "RealPos.h"
#include "RlTypedFunction.h"


namespace RevLanguage {
    
    class Func_phyloDiversity : public TypedFunction<RealPos> {
        
    public:
        Func_phyloDiversity( void );
        
        // Basic utility functions
        Func_phyloDiversity*                                                clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*                                createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif
