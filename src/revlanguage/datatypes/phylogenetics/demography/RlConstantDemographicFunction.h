#ifndef RlConstantDemographicFunction_H
#define RlConstantDemographicFunction_H

#include "ConstantDemographicFunction.h"
#include "RlDemographicFunction.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class ConstantDemographicFunction : public DemographicFunction {
        
    public:
        
        ConstantDemographicFunction(void);                                                                                                  //!< Default constructor
        
        // Basic utility functions
        virtual ConstantDemographicFunction*        clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getDemographicFunctionName(void) const;                                                 //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   theta;
        
    };
    
}

#endif
