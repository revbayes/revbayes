#ifndef Dist_beta_H
#define Dist_beta_H

#include "BetaDistribution.h"
#include "Probability.h"
#include "RlTypedDistribution.h"
#include "RlProbabilityContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the beta distribution.
     *
     * The RevLanguage wrapper of the beta distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the BetaDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-08-08, version 1.0
     *
     */
    class Dist_beta : public ProbabilityContinuousDistribution {
        
    public:
        Dist_beta( void );
        virtual ~Dist_beta();
        
        // Basic utility functions
        Dist_beta*                                      clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        virtual Probability*                            createRandomVariable(void) const;
        RevBayesCore::BetaDistribution*                 createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       beta;
        
    };
    
}

#endif
