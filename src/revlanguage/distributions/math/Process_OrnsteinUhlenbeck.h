#ifndef Process_OrnsteinUhlenbeck_H
#define Process_OrnsteinUhlenbeck_H

#include "OrnsteinUhlenbeckProcess.h"
#include "RlContinuousDistribution.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the Ornstein-Uhlenbeck process.
     *
     * The RevLanguage wrapper of the Ornstein-Uhlenbeck takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-02-23, version 1.0
     *
     */
    class OrnsteinUhlenbeckProcess : public ContinuousDistribution {
        
    public:
        OrnsteinUhlenbeckProcess( void );
        
        // Basic utility functions
        OrnsteinUhlenbeckProcess*                       clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::OrnsteinUhlenbeckProcess*         createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       x0;
        RevPtr<const RevVariable>                       theta;
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       sigma;
        RevPtr<const RevVariable>                       time;
        
    };
    
}

#endif
