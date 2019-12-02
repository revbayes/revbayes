#ifndef __revbayes_proj__Dist_cppNormal__
#define __revbayes_proj__Dist_cppNormal__

#include "CompoundPoissonNormalDistribution.h"
#include "Natural.h"
#include "RlContinuousDistribution.h"


namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Poisson distribution.
     *
     * The RevLanguage wrapper of the Poisson distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the Dist_cppNormal.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2014-03-18, version 1.0
     *
     */
    class Dist_cppNormal :  public ContinuousDistribution {
        
    public:
        Dist_cppNormal( void );                                                                                                                             //!< Default constructor
        
        // Basic utility functions
        Dist_cppNormal*                                     clone(void) const;                                                                              //!< Clone the object
        static const std::string&                           getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                        //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                        //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                                  //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::CompoundPoissonNormalDistribution*    createDistribution(void) const;                                                                 //!< Create the internal distribution object
        
    protected:
        
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                           lambda;                                                                                         //!< The rate of the distribution
        RevPtr<const RevVariable>                           mu;                                                                                             //!< The mean of the distribution
        RevPtr<const RevVariable>                           sigma;                                                                                          //!< The sd of the distribution
        
    };
    
}



#endif /* defined(__revbayes_proj__Dist_cppNormal__) */
