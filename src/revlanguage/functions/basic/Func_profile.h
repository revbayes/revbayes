/**
 * @file
 * This file contains the declaration and implementation
 * of the Func_profile, which is used to set the seed of the random number generator.
 *
 * @brief Declaration and implementation of Func_profile
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2020-03-24 10:25:58 +0200 (Tue, 12 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_profile.h 1626 2020-03-24 08:25:58Z meyerx $
 */

#ifndef FUNC_PROFILE_H
#define FUNC_PROFILE_H

#include <iosfwd>
#include <vector>

#include "Procedure.h"
#include "RevPtr.h"

namespace RevBayesCore { class RbHelpReference; }

namespace RevLanguage {
class ArgumentRules;
class RevVariable;
class TypeSpec;

    class Func_profile : public Procedure {

    public:
        Func_profile();

        // Basic utility functions
        Func_profile*                                   clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get language type of the object

        // Regular functions
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                  //!< Get type of return value


        RevPtr<RevVariable>                             execute(void);                                              //!< Execute function

    protected:

    };

}


#endif
