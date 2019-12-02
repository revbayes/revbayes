#ifndef Dist_bimodalNorm_H
#define Dist_bimodalNorm_H

#include "BimodalNormalDistribution.h"
#include "RlContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the bimodal normal distribution.
     *
     * The RevLanguage wrapper of the bimodal normal distribution takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-04-09, version 1.0
     *
     */
    class Dist_bimodalNorm :  public ContinuousDistribution {
        
    public:
        Dist_bimodalNorm( void );
        
        // Basic utility functions
        Dist_bimodalNorm*                               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::BimodalNormalDistribution*        createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       mean1;
        RevPtr<const RevVariable>                       mean2;
        RevPtr<const RevVariable>                       sd1;
        RevPtr<const RevVariable>                       sd2;
        RevPtr<const RevVariable>                       prob;
        
    };
    
}

#endif
