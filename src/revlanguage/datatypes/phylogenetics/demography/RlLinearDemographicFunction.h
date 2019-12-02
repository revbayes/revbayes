#ifndef RlLinearDemographicFunction_H
#define RlLinearDemographicFunction_H

#include "LinearDemographicFunction.h"
#include "RlDemographicFunction.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class LinearDemographicFunction : public DemographicFunction {
        
    public:
        
        LinearDemographicFunction(void);                                                                                                  //!< Default constructor
        
        // Basic utility functions
        virtual LinearDemographicFunction*          clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getDemographicFunctionName(void) const;                                                 //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   theta_recent;
        RevPtr<const RevVariable>                   theta_ancient;
        RevPtr<const RevVariable>                   time_recent;
        RevPtr<const RevVariable>                   time_ancient;
        
    };
    
}

#endif
