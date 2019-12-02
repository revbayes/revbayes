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
    
        std::vector<std::string>                        getHelpAuthor(void) const;                                                              //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                                                         //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                                             //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                                             //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                                                          //!< Get some references/citations for this function
        std::vector<std::string>                        getHelpSeeAlso(void) const;                                                             //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                                               //!< Get the title of this help entry

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    
    private:
        
        RevPtr<const RevVariable>                       mean;
        RevPtr<const RevVariable>                       covariance;
        RevPtr<const RevVariable>                       precision;
        RevPtr<const RevVariable>                       scale;
        
    };
    
}

#endif
