#ifndef Func_coala_H
#define Func_coala_H

#include "ModelVector.h"
#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the coala function.
     *
     * The RevLanguage wrapper of the coala function connects
     * the variables/parameters of the function and creates the internal coalaFunction object.
     * Please read the coalaFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_coala : public TypedFunction< ModelVector<Real> > {
        
    public:
        Func_coala();
        
        // Basic utility functions
        Func_coala*                                                                 clone(void) const;                                          //!< Clone the object
        static const std::string&                                                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::RbVector< double > >*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                        getArgumentRules(void) const;                               //!< Get argument rules
        
        
    protected:
        
        std::vector<std::string>                                                    getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                                                 getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                                                 getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                                                 getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>                                  getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                                                    getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                                                 getHelpTitle(void) const;                                   //!< Get the title of this help entry

    };
    
}


#endif


