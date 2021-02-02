#ifndef Dist_Cauchy_H
#define Dist_Cauchy_H

#include "CauchyDistribution.h"
#include "RlContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the half-Cauchy distribution.
     *
     * The RevLanguage wrapper of the half-Cauchy distribution takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author Andrew Magee
     * @since 2017-02-04, version 1.0
     *
     */
    class Dist_Cauchy : public ContinuousDistribution {
    
    public:
        Dist_Cauchy( void );
    
        // Basic utility functions
        Dist_Cauchy*                                clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                                     //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                        //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                        //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                                  //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')

    
        // Distribution functions you have to override
        RevBayesCore::CauchyDistribution*               createDistribution(void) const;
    
    protected:
    
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    
    private:
        
        RevPtr<const RevVariable>                          location;
        RevPtr<const RevVariable>                          scale;
    
    };
    
}

#endif
