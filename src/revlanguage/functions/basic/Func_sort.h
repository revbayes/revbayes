/**
 * @file
 * This file contains the declaration and implementation
 * of the Func_sort, which is used to set the seed of the random number generator.
 *
 * @brief Declaration and implementation of Func_sort
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-12 10:25:58 +0200 (Tue, 12 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_sort.h 1626 2012-06-12 08:25:58Z hoehna $
 */

#ifndef Func_sort_H
#define Func_sort_H

#include "ModelVector.h"
#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    class Func_sort : public TypedFunction< ModelVector< Real > > {
        
    public:
        Func_sort();
        
        // Basic utility functions
        Func_sort*                  clone(void) const;                                          //!< Clone the object
        static const std::string&   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&             getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* createFunction( void ) const;
        const ArgumentRules&        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif


