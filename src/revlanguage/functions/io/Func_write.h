/**
 * @file
 * This file contains the declaration of Func_write, which writes the contents of the variables into a file.
 *
 * @brief Declaration of Func_write
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-27 16:10:36 +0200 (Sun, 27 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_readDiscreteCharacterData.h 1587 2012-05-27 14:10:36Z hoehna $
 */

#ifndef Func_write_H
#define Func_write_H

#include "Procedure.h"


namespace RevLanguage {
    
    class Func_write : public Procedure {
        
    public:
        Func_write();
        
        // Basic utility functions
        Func_write*                             clone(void) const;                                          //!< Clone the object
        static const std::string&               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::vector<std::string>                getFunctionNameAliases(void) const;                         //!< Get aliases for the Rev name of the function
        std::string                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                         getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                     execute(void);                                              //!< Execute function
        const ArgumentRules&                    getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                         getReturnType(void) const;                                  //!< Get type of return value
        
    private:
        
        int                                     process_ID;
    };
    
}

#endif

