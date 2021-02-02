/**
 * @file
 * This file contains the declaration of the RevLanguage wrapper of the AncestralStateMonitor class.
 *
 * @brief Declaration of Mntr_AncestralState
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Mntr_AncestralState.h $
 */

#ifndef Mntr_AncestralState_H
#define Mntr_AncestralState_H

#include <ostream>
#include <vector>

#include "RlMonitor.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
class TypeSpec;
    
    class Mntr_AncestralState : public Monitor {
        
    public:
        
        Mntr_AncestralState(void);                                                                                                          //!< Default constructor
        
        // Basic utility functions
        virtual Mntr_AncestralState*                clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMonitorName(void) const;                                                             //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
		std::vector<RevPtr<const RevVariable> >     vars;
        RevPtr<const RevVariable>                   filename;
        RevPtr<const RevVariable>                   printgen;
		RevPtr<const RevVariable>                   tree;
		RevPtr<const RevVariable>                   ctmc;
        RevPtr<const RevVariable>                   separator;
        RevPtr<const RevVariable>                   append;
        RevPtr<const RevVariable>                   version;
		RevPtr<const RevVariable>                   monitorType;
    };
    
}


#endif
