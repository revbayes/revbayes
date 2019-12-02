//
//  Dist_varianceGamma.h
//  revbayes-proj
//
//  Created by Michael Landis on 2/8/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Dist_varianceGamma_h
#define Dist_varianceGamma_h

#include "VarianceGammaDistribution.h"
#include "RlContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the normal distribution.
     *
     * The RevLanguage wrapper of the normal distribution takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-08-06, version 1.0
     *
     */
    class Dist_varianceGamma : public ContinuousDistribution {
        
    public:
        Dist_varianceGamma( void );
        
        // Basic utility functions
        Dist_varianceGamma*                             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                                     //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                        //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                        //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                                  //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::VarianceGammaDistribution*        createDistribution(void) const;
        
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
        
        RevPtr<const RevVariable>                          mu;
        RevPtr<const RevVariable>                          kappa;
        RevPtr<const RevVariable>                          tau;
        RevPtr<const RevVariable>                          time;
        
    };
    
}

#endif /* Dist_varianceGamma_h */
