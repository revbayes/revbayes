/**
 * @file
 * This file contains the declaration of Func_ls, which is
 * the function used to list the content of the workspace.
 *
 * @brief Declaration of Func_ls
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

#ifndef Func_ls_H
#define Func_ls_H

#include "Procedure.h"

namespace RevLanguage {
    
//    class TypeTable;
    class Func_ls : public Procedure {
        
    public:
        Func_ls( void );
        
        // Basic utility functions
        Func_ls*                                        clone(void) const;                                                          //!< Clone object
        static const std::string&                       getClassType(void);                                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                    //!< Get language type of the object
        
        // Func_ls functions
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

    private:
//        std::map<std::string, RevObject*>               applyTypeTableFilters( const std::map<std::string, RevObject*>& table, std::string filter ) const;
//       const std::multimap<std::string, Function*>
        std::string                                     makeFilteredFunctionString( const std::multimap<std::string, Function*>& table, std::string filter ) const;
        void                                            printTypes(bool printAll, std::string filter) const;
        void                                            printVariables(bool printAll, std::string filter) const;
        void                                            printFunctions(bool printAll, std::string filter) const;
        bool                                            isTokenCategorized( std::string name ) const;
        bool                                            matchToken( std::string name, std::string filter ) const;
        bool                                            matchAnalysisToken( std::string name ) const;
        bool                                            matchMathToken( std::string name ) const;
        bool                                            matchFileToken( std::string name ) const;
        bool                                            matchSummaryToken( std::string name ) const;
        bool                                            matchUtilityToken( std::string name ) const;
    };
    
}

#endif

