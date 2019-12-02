//
//  Mntr_CharacterHistorySummary.hpp
//  revbayes-proj
//
//  Created by Michael Landis on 12/3/18.
//  Copyright © 2018 Michael Landis. All rights reserved.
//

#ifndef Mntr_CharacterHistorySummary_h
#define Mntr_CharacterHistorySummary_h

#include "CharacterHistorySummaryMonitor.h"
#include "RlMonitor.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Mntr_CharacterHistorySummary : public Monitor {
        
    public:
        
        Mntr_CharacterHistorySummary(void);                                                                                              //!< Default constructor (0.0)
        
        // Basic utility functions
        virtual Mntr_CharacterHistorySummary*    clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal Mntr_CharacterHistorySummary.
        static const std::string&                   getClassType(void);                                                                     //!< Get class name
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMonitorName(void) const;                                                             //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        std::set<RevPtr<const RevVariable> >        vars;
        RevPtr<const RevVariable>                   filename;
        RevPtr<const RevVariable>                   tree;
        RevPtr<const RevVariable>                   ctmc;
        RevPtr<const RevVariable>                   printgen;
        RevPtr<const RevVariable>                   separator;
        RevPtr<const RevVariable>                   prior;
        RevPtr<const RevVariable>                   posterior;
        RevPtr<const RevVariable>                   likelihood;
        RevPtr<const RevVariable>                   type;
        //        RevPtr<const RevVariable>                      counts;
        //        RevPtr<const RevVariable>                      events;
        RevPtr<const RevVariable>                      style;
        RevPtr<const RevVariable>                      append;
        
    };
    
}

#endif /* Mntr_CharacterHistorySummary_h */
