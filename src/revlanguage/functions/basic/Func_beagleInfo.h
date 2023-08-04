/**
 * @file
 * This file contains the declaration of Func_beagleInfo, which is
 * used to provide the user with information on the available
 * BEAGLE resources and options.
 *
 * @brief Declaration of Func_beagleInfo
 *
 * @author Daniel L. Ayres
 * @license GPL version 3
 */

#ifndef Func_beagleInfo_h
#define Func_beagleInfo_h

#include "Procedure.h"

namespace RevLanguage {
    
    class Func_beagleInfo : public Procedure {
        
    public:
        Func_beagleInfo( void );
        
        // Basic utility functions
        Func_beagleInfo*                              clone(void) const;                 //!< Clone object
        static const std::string&                     getClassType(void);                //!< Get Rev type
        static const TypeSpec&                        getClassTypeSpec(void);            //!< Get class type spec
        std::string                                   getFunctionName(void) const;
        const TypeSpec&                               getTypeSpec(void) const;           //!< Get language type of the object
        
        // Func_beagleInfo functions
        const ArgumentRules&                          getArgumentRules(void) const;      //!< Get argument rules
        const TypeSpec&                               getReturnType(void) const;         //!< Get type of return val
        
        RevPtr<RevVariable>                           execute(void);                     //!< Execute function
        
    protected:
        
        std::vector<std::string>                      getHelpAuthor(void) const;         //!< Get the author(s) of this function
        std::string                                   getHelpDescription(void) const;    //!< Get the description for this function
        std::string                                   getHelpDetails(void) const;        //!< Get the more detailed description of the function
        std::string                                   getHelpExample(void) const;        //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>    getHelpReferences(void) const;     //!< Get some references/citations for this function
        std::vector<std::string>                      getHelpSeeAlso(void) const;        //!< Get suggested other functions
        std::string                                   getHelpTitle(void) const;          //!< Get the title of this help entry

    };
    
}

#endif /* Func_beagleInfo_h */
