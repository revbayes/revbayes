//
//  Mntr_JointConditionalJointConditionalAncestralState.h
//  revbayes-proj
//
//  Created by Michael Landis on 2/28/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef Mntr_JointConditionalJointConditionalAncestralState_H
#define Mntr_JointConditionalJointConditionalAncestralState_H


#include <ostream>
#include <vector>

#include "RlMonitor.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
class TypeSpec;
    
    class Mntr_JointConditionalAncestralState : public Monitor {
        
    public:
        
        Mntr_JointConditionalAncestralState(void);                                                                                              //!< Default constructor
        
        // Basic utility functions
        virtual Mntr_JointConditionalAncestralState*    clone(void) const;                                                                      //!< Clone object
        void                                            constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getMonitorName(void) const;                                                             //!< Get the name used for the constructor function in Rev.
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                         getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                    printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
		std::vector<RevPtr<const RevVariable> >         vars;
        RevPtr<const RevVariable>                       filename;
        RevPtr<const RevVariable>                       printgen;
        RevPtr<const RevVariable>                       version;
		RevPtr<const RevVariable>                       tree;
		RevPtr<const RevVariable>                       ctmc;
        RevPtr<const RevVariable>                       cdbdp;
        RevPtr<const RevVariable>                       glhbdsp;
        RevPtr<const RevVariable>                       separator;
        RevPtr<const RevVariable>                       append;
		RevPtr<const RevVariable>                       monitorType;
        RevPtr<const RevVariable>                       withTips;
        RevPtr<const RevVariable>                       withStartStates;
    };
    
}

#endif
