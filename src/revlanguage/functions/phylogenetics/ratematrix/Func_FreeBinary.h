//
//  Func_FreeBinary.h
//  rb_mlandis
//
//  Created by Michael Landis on 9/19/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

/**
 * @file
 * This file contains the declaration of the RevLanguage Jukes-Cantor (JC) function, which
 * is used to create deterministic variable associated with the jc function.
 *
 * @brief Declaration and implementation of RlFunc_FreeBinary
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func__add.h 1406 2012-04-20 02:06:14Z hoehna $
 */


#ifndef Func_FreeBinary_H
#define Func_FreeBinary_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_FreeBinary : public TypedFunction<RateMatrix> {
        
    public:
        Func_FreeBinary( void );
        
        // Basic utility functions
        Func_FreeBinary*                                                    clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >*         createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
