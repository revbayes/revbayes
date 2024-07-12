#ifndef Dist_Inverse_H
#define Dist_Inverse_H

// Need to include headers for all possible values of valType
#include "Integer.h"
#include "Natural.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "Simplex.h"
#include "IidDistribution.h"
#include "ModelVector.h"
#include "TypedDistribution.h"
#include "InverseDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {

    template<typename valType>
    class Dist_Inverse : public TypedDistribution< valType > {
        
    public:
        Dist_Inverse( void );
        virtual ~Dist_Inverse();
        
        // Basic utility functions
        Dist_Inverse*                                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::InverseDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       dist;
    };
    
}

#endif // Dist_Inverse_H
