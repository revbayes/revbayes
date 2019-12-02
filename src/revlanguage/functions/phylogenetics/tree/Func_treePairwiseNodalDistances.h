/**
 * @file
 * This file contains the declaration of the RevLanguage pairwise distances function, which
 * is used to create a deterministic variable associated with the pairwise nodal distances matrix.
 *
 * @brief Declaration and implementation of Func_treePairwiseNodalDistances
 *
 *
 *  Func_treePairwiseNodalDistances.h Will Freyman 5/21/2018
 */


#ifndef Func_treePairwiseNodalDistances_H
#define Func_treePairwiseNodalDistances_H

#include "RlDistanceMatrix.h"
#include "RlTypedFunction.h"

#include <string>


namespace RevLanguage {

    class Func_treePairwiseNodalDistances : public TypedFunction< DistanceMatrix > {
        
    public:
        Func_treePairwiseNodalDistances( void );
        
        // Basic utility functions
        Func_treePairwiseNodalDistances*                  					clone(void) const;                                          //!< Clone the object
        static const std::string&                       					getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          					getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 					getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::DistanceMatrix >*        createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            					getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
