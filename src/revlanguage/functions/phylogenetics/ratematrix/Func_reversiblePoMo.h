/**
 * @file
 * This file contains the declaration of the RevLanguage pomo function, which
 * is used to created a deterministic variable associated with the pomo function.
 *
 * @brief Declaration and implementation of Func_reversiblePoMo
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_reversiblePoMo.h 1406 2014-08-22 02:06:14Z Boussau $
 */


#ifndef Func_reversiblePoMo_H
#define Func_reversiblePoMo_H

#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {

    class Func_reversiblePoMo : public TypedFunction<RateMatrix> {

    public:
        Func_reversiblePoMo( void );

        // Basic utility functions
        Func_reversiblePoMo*                                            clone(void) const;                                          //!< Clone the object
        static const std::string&                                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RateGenerator>*      createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                               //!< Get argument rules

    };

}

#endif
