#ifndef Dist_gamma_H
#define Dist_gamma_H

#include "GammaDistribution.h"
#include "RlPositiveContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the gamma distribution.
     *
     * The RevLanguage wrapper of the gamma distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the GammaDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-08-08, version 1.0
     *
     */
    class Dist_gamma : public PositiveContinuousDistribution {
        
    public:
        Dist_gamma( void );
        virtual ~Dist_gamma();
        
        // Basic utility functions
        Dist_gamma*                                     clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::GammaDistribution*                createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       shape;
        RevPtr<const RevVariable>                       rate;
        
    };
    
}

#endif
