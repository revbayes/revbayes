#ifndef RlMultiValueEvent_H
#define RlMultiValueEvent_H

#include <string>
#include <vector>
#include <iosfwd>

#include "ModelObject.h"
#include "MultiValueEvent.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class RbHelpReference; }


namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;
    
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
        
        void                                        initMethods(void);
        //        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
        
        
    };
}

#endif


