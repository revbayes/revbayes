#ifndef Dist_Laplaceal_H
#define Dist_Laplaceal_H

#include "LaplaceDistribution.h"
#include "RlContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the Laplace distribution.
     *
     * The RevLanguage wrapper of the Laplace distribution takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     n
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team 
     * @since 2017-10-19
     *
     */
    class Dist_Laplace : public ContinuousDistribution {
    
    public:
        Dist_Laplace( void );
    
        Dist_Laplace*                                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        RevBayesCore::LaplaceDistribution*              createDistribution(void) const;
    
    protected:
    
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    private:
        
        RevPtr<const RevVariable>                       mean;
        RevPtr<const RevVariable>                       scale;
    
    };
    
}

#endif
