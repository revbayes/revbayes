#ifndef Func_writeCountsFile_H
#define Func_writeCountsFile_H

#include "Procedure.h"


namespace RevLanguage {

    /**
     * Function that takes in character data and writes it into a file in CountsFile format.
     *
     *
     * This function is a very basic function that simply takes in a character data object
     * and delegate the writing/outputting of the data into a file to the CountsFileWrite.
     * See CountsFileWriter.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class Func_writeCountsFile : public Procedure {
        
    public:
        // Basic utility functions
        Func_writeCountsFile*               clone(void) const;                                          //!< Clone the object
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

#endif

