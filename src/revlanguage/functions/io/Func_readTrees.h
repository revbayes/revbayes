#ifndef Func_readTrees_H
#define Func_readTrees_H

#include "Procedure.h"
#include "RbFileManager.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {

    class Func_readTrees : public Procedure {

    public:
        // Basic utility functions
        Func_readTrees*                                         clone(void) const;                                          //!< Clone the object
        static const std::string&                               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                         getTypeSpec(void) const;                                    //!< Get language type of the object

        // Regular functions
        RevPtr<RevVariable>                                     execute(void);                                              //!< Execute function
        const ArgumentRules&                                    getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                         getReturnType(void) const;                                  //!< Get type of return value

    protected:

        std::vector<std::string>                                getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                             getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                             getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                             getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>              getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                                getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                             getHelpTitle(void) const;                                   //!< Get the title of this help entry



    };

}

#endif
