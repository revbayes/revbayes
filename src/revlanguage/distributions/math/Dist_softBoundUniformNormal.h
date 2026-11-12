#ifndef Dist_SoftBoundUniformNormal_H
#define Dist_SoftBoundUniformNormal_H

#include "SoftBoundUniformNormalDistribution.h"
#include "RlContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the uniform distribution with normal soft bounds.
     *
     * The RevLanguage wrapper of the uniform distribution with normal soft bounds takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-20, version 1.0
     *
     */
    class Dist_SoftBoundUniformNormal : public ContinuousDistribution {
        
    public:
        Dist_SoftBoundUniformNormal( void );
        
        // Basic utility functions
        Dist_SoftBoundUniformNormal*                    clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::SoftBoundUniformNormalDistribution*   createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       min;
        RevPtr<const RevVariable>                       max;
        RevPtr<const RevVariable>                       sd;
        RevPtr<const RevVariable>                       prob;
        RevPtr<const RevVariable>                       boundary;

    };
    
}

#endif
