//
//  Func_formatDiscreteCharacterData.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/1/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Func_formatDiscreteCharacterData_h
#define Func_formatDiscreteCharacterData_h


#include "Procedure.h"


namespace RevLanguage {
    
    /**
     * Function that writes character data as a delimited text file.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Will Freyman)
     * @since 2013-04-15, version 1.0
     */
    class Func_formatDiscreteCharacterData : public Procedure {
        
    public:
        // Basic utility functions
        Func_formatDiscreteCharacterData*                    clone(void) const;                                          //!< Clone the object
        static const std::string&           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                 execute(void);                                              //!< Execute function
        const ArgumentRules&                getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                     getReturnType(void) const;                                  //!< Get type of return value
        
        
    };
    
}

#endif /* Func_formatDiscreteCharacterData_h */
