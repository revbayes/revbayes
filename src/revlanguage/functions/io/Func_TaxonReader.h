#ifndef Func_TaxonReader_H
#define Func_TaxonReader_H

#include "Procedure.h"

#include <string>


namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the TaxonReader function.
     *
     * The RevLanguage wrapper of the TaxonReader function which reads in a set of taxa specification from a delimitted file.
     * Please read the TaxonReader.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-01, version 1.0
     *
     */
    class Func_TaxonReader :  public Procedure {
        
    public:
        // Basic utility functions
        Func_TaxonReader*                   clone(void) const;                                          //!< Clone the object
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

