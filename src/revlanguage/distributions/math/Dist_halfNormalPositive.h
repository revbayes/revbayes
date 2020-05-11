#ifndef Dist_halfNormalPositive_H
#define Dist_halfNormalPositive_H

#include "HalfNormalDistribution.h"
#include "RlPositiveContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the half-normal distribution.
     *
     * The RevLanguage wrapper of the half-normal distribution takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author Andrew Magee
     * @since 2017-03-29, version 1.0
     *
     */
    class Dist_halfNormalPositive :  public PositiveContinuousDistribution {
        
    public:
        Dist_halfNormalPositive( void );
        
        // Basic utility functions
        Dist_halfNormalPositive*                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::HalfNormalDistribution*      createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       offset;
        RevPtr<const RevVariable>                       sd;
        
    };
    
}

#endif
