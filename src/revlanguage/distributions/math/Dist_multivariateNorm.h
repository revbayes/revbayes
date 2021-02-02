#ifndef Dist_multivariateNorm_H
#define Dist_multivariateNorm_H

#include "MultivariateNormalDistribution.h"
#include "RlTypedDistribution.h"
#include "ModelVector.h"
#include "Real.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the multivariate normal distribution.
     *
     * The RevLanguage wrapper of the multivariate normal distribution takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John Huelsenbeck and Risa Takenaka)
     * @since 2012-08-06, version 1.0
     *
     */
    class Dist_multivariateNorm :  public TypedDistribution<ModelVector<Real> > {
    
    public:
                                                        Dist_multivariateNorm( void );
    
        // Basic utility functions
        Dist_multivariateNorm*                          clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
//        // Member method functions
//        virtual RevPtr<RevVariable>                     executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
//        virtual const MethodTable&                      getMethods(void) const;                                                                     //!< Get member methods

    
        // Distribution functions you have to override
        RevBayesCore::MultivariateNormalDistribution*   createDistribution(void) const;
    
    protected:
    
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    
    private:
        
        RevPtr<const RevVariable>                       mean;
        RevPtr<const RevVariable>                       covariance;
        RevPtr<const RevVariable>                       precision;
        RevPtr<const RevVariable>                       scale;
        
    };
    
}

#endif
