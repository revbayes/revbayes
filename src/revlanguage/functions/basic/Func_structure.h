#ifndef Func_structure_H
#define Func_structure_H

#include "Procedure.h"

namespace RevLanguage {
    
    /**
     * @brief Rev function to print the structural information of a workspace variable.
     * 
     * The structure function is used by Rev to print structural information,
     * such as children and parents as well type and value of a variable.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-01-21
     *
     */
    class Func_structure : public Procedure {
        
    public:
        Func_structure( void );                                                                             //!< Default constructor
        
        // Basic utility functions
        Func_structure*                                 clone(void) const;                                          //!< Clone object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::vector<std::string>                        getFunctionNameAliases(void) const;                         //!< Get aliases for the Rev name of the function
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Func_source functions
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                  //!< Get type of return val
        
        RevPtr<RevVariable>                             execute(void);                                              //!< Execute function
        
    protected:
        
        std::vector<std::string>                        getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                        getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                   //!< Get the title of this help entry

    };
    
}

#endif

