#ifndef Mntr_StochasticVariable_H
#define Mntr_StochasticVariable_H

#include "StochasticVariableMonitor.h"
#include "RlFileMonitor.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Mntr_StochasticVariable : public FileMonitor {
        
    public:
        
        Mntr_StochasticVariable(void);                                                                                                  //!< Default constructor (0.0)
        
        // Basic utility functions
        virtual Mntr_StochasticVariable*            clone(void) const;                                                                  //!< Clone object
        void                                        constructInternalObject(void);                                                      //!< We construct the a new internal Mntr_StochasticVariable.
        static const std::string&                   getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                             //!< Get class type spec
        std::string                                 getMonitorName(void) const;                                                         //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                      //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                            //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                  //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable

    };
    
}

#endif
