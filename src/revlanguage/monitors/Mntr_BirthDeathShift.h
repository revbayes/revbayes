#ifndef Mntr_BirthDeathShift_H
#define Mntr_BirthDeathShift_H

#include <ostream>
#include <vector>

#include "RlFileMonitor.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"


namespace RevLanguage {
class TypeSpec;
    
    class Mntr_BirthDeathShift : public FileMonitor {
        
    public:
        
        Mntr_BirthDeathShift(void);                                                                                            //!< Default constructor
        
        // Basic utility functions
        virtual Mntr_BirthDeathShift*              clone(void) const;                                                                      //!< Clone object
        void                                            constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getMonitorName(void) const;                                                             //!< Get the name used for the constructor function in Rev.
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                         getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                    printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                       bdsp;

    };
    
}


#endif
