#ifndef Dist_Log_H
#define Dist_Log_H

#include "IidDistribution.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TransformedDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    class Dist_Log : public TypedDistribution<RealPos> {
        
    public:
        Dist_Log( void );
        virtual ~Dist_Log();
        
        // Basic utility functions
        Dist_Log*                                       clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::TransformedDistribution*          createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       log_distribution;
        
    };
    
}

#endif // Dist_Log_H
