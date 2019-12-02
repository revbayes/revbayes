#ifndef Mntr_Probability_H
#define Mntr_Probability_H

#include "ModelMonitor.h"
#include "RlMonitor.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Mntr_Probability : public Monitor {
        
    public:
        
        Mntr_Probability(void);                                                                                                                   //!< Default constructor
        
        // Basic utility functions
        virtual Mntr_Probability*                   clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMonitorName(void) const;                                                             //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   filename;
        RevPtr<const RevVariable>                   printgen;
        RevPtr<const RevVariable>                   separator;
        RevPtr<const RevVariable>                   prior;
        RevPtr<const RevVariable>                   posterior;
        RevPtr<const RevVariable>                   likelihood;
        RevPtr<const RevVariable>                   append;
        RevPtr<const RevVariable>                   version;
        
    };
    
}

#endif
