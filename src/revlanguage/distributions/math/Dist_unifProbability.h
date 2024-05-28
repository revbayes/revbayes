#ifndef Dist_unifProbability_H
#define Dist_unifProbability_H

#include "UniformDistribution.h"
#include "RlProbabilityContinuousDistribution.h"
#include "RlTypedDistribution.h"
#include "Probability.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the uniform distribution on probabilities.
     *
     * The RevLanguage wrapper of the uniform distribution on probabilities simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * We need this second implementation because the return type must be able to guarantee to be a probability
     * if both boundaries are probabilities.
     * See the UniformDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-15, version 1.0
     *
     */
    class Dist_unifProbability : public ProbabilityContinuousDistribution {
        
    public:
        Dist_unifProbability( void );
        virtual ~Dist_unifProbability();
        
        // Basic utility functions
        Dist_unifProbability*                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::UniformDistribution*              createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       lower;
        RevPtr<const RevVariable>                       upper;
        
    };
    
}

#endif
