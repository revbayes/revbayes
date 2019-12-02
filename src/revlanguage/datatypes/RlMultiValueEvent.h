#ifndef RlMultiValueEvent_H
#define RlMultiValueEvent_H

#include "ModelObject.h"
#include "MultiValueEvent.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of a MultiValueEvent, which is
     * the class that contains the name of a MultiValueEvent with associated species and date.
     *
     * @brief Declaration of MultiValueEvent
     *
     * (c) Copyright 2009-
     * @date Last modified: $Date: $
     * @author The RevBayes Development Core Team
     * @license GPL version 3
     *
     *
     */
    
    class MultiValueEvent : public ModelObject<RevBayesCore::MultiValueEvent> {
        
    public:
        MultiValueEvent(void);                                                                        //!< Constructor requires character type
        MultiValueEvent(RevBayesCore::MultiValueEvent *v);                                                      //!< Constructor requires character type
        MultiValueEvent(const RevBayesCore::MultiValueEvent &v);                                                //!< Constructor requires character type
        MultiValueEvent(RevBayesCore::TypedDagNode<RevBayesCore::MultiValueEvent> *n);                          //!< Constructor requires character type
        
        typedef RevBayesCore::MultiValueEvent         valueType;
        
        // Basic utility functions
        MultiValueEvent*                            clone(void) const;                                                                  //!< Clone object
        //        void                                        constructInternalObject(void);                                                      //!< We construct the a new internal MCMC object.
        static const std::string&                   getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                             //!< Get class type spec
//        const MemberRules&                          getParameterRules(void) const;                                                      //!< Get member rules (const)
        const TypeSpec&                             getTypeSpec(void) const;                                                            //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f); //!< Map member methods to internal functions
        
        std::string                                 getGuiName(void) { return "MultiValueEvent"; }
        std::string                                 getGuiUnicodeSymbol(void) { return "OTU"; }
        std::string                                 getGuiInfo(void) { return ""; }
        
    protected:
        
        std::vector<std::string>                    getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                 getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                 getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                 getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>  getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                    getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                 getHelpTitle(void) const;                                   //!< Get the title of this help entry
        
        void                                        initMethods(void);
        //        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
        
        
    };
}

#endif


