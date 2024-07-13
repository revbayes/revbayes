#ifndef Transform_Exp_H
#define Transform_Exp_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "LogDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    class Transform_Exp : public TypedDistribution<RealPos> {
        
    public:
        Transform_Exp( void );
        virtual ~Transform_Exp();
        
        // Basic utility functions
        Transform_Exp*                                  clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::LogDistribution*                  createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:
        RevPtr<const RevVariable>                       base_distribution;
        
    };
    
}

#endif // Transform_Exp_H
