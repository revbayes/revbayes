#ifndef Dist_inverseGamma_H
#define Dist_inverseGamma_H

#include "InverseGammaDistribution.h"
#include "RlPositiveContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the inverse Gamma distribution.
     *
     * The RevLanguage wrapper of the inverse Gamma distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the inverseGammaDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2017-03-01, version 1.0
     *
     */
    class Dist_inverseGamma :  public PositiveContinuousDistribution {
        
    public:
        Dist_inverseGamma( void );
        virtual ~Dist_inverseGamma();
        
        // Basic utility functions
        Dist_inverseGamma*                              clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::InverseGammaDistribution*         createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       shape;
        RevPtr<const RevVariable>                       rate;
        
    };
    
}

#endif
