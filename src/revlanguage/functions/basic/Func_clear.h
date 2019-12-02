/**
 * @file
 * This file contains the declaration of Func_clear, which is
 * the function used to clear the workspace.
 *
 * @brief Declaration of Func_clear
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-22 12:06:58 -0700 (Sun, 22 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @interface RbFunction
 * @package functions
 * @since Version 1.0, 2012-09-07
 *
 * $Id: Func_source.h 1431 2012-04-22 19:06:58Z hoehna $
 */

#ifndef Func_clear_H
#define Func_clear_H

#include "Procedure.h"

namespace RevLanguage {
    
    class Func_clear : public Procedure {
        
    public:
        Func_clear( void );
        
        // Basic utility functions
        Func_clear*                                     clone(void) const;                                                          //!< Clone object
        static const std::string&                       getClassType(void);                                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                    //!< Get language type of the object
        
        // Func_clear functions
        const ArgumentRules&                            getArgumentRules(void) const;                                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                                  //!< Get type of return val
        
        RevPtr<RevVariable>                             execute(void);                                                              //!< Execute function
        
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

