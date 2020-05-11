/**
 * @file
 * This file contains the declaration of Func__range, which is used
 * to look up a variable.
 *
 * @brief Declaration of Func__range
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-10 13:42:43 +0200 (Thu, 10 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func__range.h 1518 2012-05-10 11:42:43Z hoehna $
 */

#ifndef Func_range_H
#define Func_range_H

#include "Procedure.h"

#include <map>
#include <string>
#include <vector>

namespace RevLanguage {

    class Func_range : public Procedure {
    
    public:
        Func_range( void );
    
        // Basic utility functions
        Func_range*                                     clone(void) const;                                          //!< Clone object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get language type of the object
    
        // Func_range functions
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                  //!< Get type of return val
    
        RevPtr<RevVariable>                             execute(void);                                              //!< Execute function
        
    protected:
        
    };
    
}

#endif

