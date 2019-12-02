#ifndef Dist_halfCauchyPositive_H
#define Dist_halfCauchyPositive_H

#include "HalfCauchyDistribution.h"
#include "RlPositiveContinuousDistribution.h"

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
    class Dist_halfCauchyPositive : public PositiveContinuousDistribution {
    
    public:
        Dist_halfCauchyPositive( void );
    
        // Basic utility functions
        Dist_halfCauchyPositive*                                clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                                     //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                        //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                        //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                                  //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')

    
        // Distribution functions you have to override
        RevBayesCore::HalfCauchyDistribution*               createDistribution(void) const;
    
    protected:
    
        std::vector<std::string>                        getHelpAuthor(void) const;                                                              //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                                                         //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                                             //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                                             //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                                                          //!< Get some references/citations for this function
//        std::vector<std::string>                        getHelpSeeAlso(void) const;                                                             //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                                               //!< Get the title of this help entry

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    
    private:
        
        RevPtr<const RevVariable>                          location;
        RevPtr<const RevVariable>                          scale;
    
    };
    
}

#endif
