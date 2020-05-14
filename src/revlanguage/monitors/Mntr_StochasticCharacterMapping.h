//
//  Mntr_StochasticCharacterMapping.h
//  RevBayes_development_branch
//
//  Created by will freyman on 12/13/16.
//  Copyright (c) 2016 will freyman. All rights reserved.
//

#ifndef __RevBayes_development_branch__Mntr_StochasticCharacterMapping__
#define __RevBayes_development_branch__Mntr_StochasticCharacterMapping__

#include <ostream>
#include <vector>

#include "RlFileMonitor.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"


namespace RevLanguage {
class TypeSpec;

    class Mntr_StochasticCharacterMapping : public FileMonitor {

    public:

        Mntr_StochasticCharacterMapping(void);                                                                                            //!< Default constructor

        // Basic utility functions
        virtual Mntr_StochasticCharacterMapping*        clone(void) const;                                                                      //!< Clone object
        void                                            constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getMonitorName(void) const;                                                             //!< Get the name used for the constructor function in Rev.
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                         getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                    printValue(std::ostream& o) const;                                                      //!< Print value (for user)

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

        RevPtr<const RevVariable>                       cdbdp;
        RevPtr<const RevVariable>                       glhbdsp;
        RevPtr<const RevVariable>                       ctmc;
        RevPtr<const RevVariable>                       include_simmap;
        RevPtr<const RevVariable>                       use_simmap_default;
        RevPtr<const RevVariable>                       index;

    };

}


#endif
