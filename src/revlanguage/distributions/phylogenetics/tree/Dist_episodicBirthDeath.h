#ifndef Dist_episodicBirthDeath_H
#define Dist_episodicBirthDeath_H

#include "EpisodicBirthDeathProcess.h"
#include "RlBirthDeathProcess.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the episodic birth-death process.
     *
     * The RevLanguage wrapper of the episodic birth-death process connects
     * the variables/parameters of the process and creates the internal Dist_episodicBirthDeath object.
     * Please read the Dist_episodicBirthDeath.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2016-02-08, version 1.0
     *
     */
    class Dist_episodicBirthDeath : public BirthDeathProcess {
        
    public:
        Dist_episodicBirthDeath( void );
        
        // Basic utility functions
        Dist_episodicBirthDeath*                        clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::EpisodicBirthDeathProcess*        createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       lambda_rates;                                                                           //!< The speciation rates
        RevPtr<const RevVariable>                       lambda_times;                                                                           //!< The speciation rate change times
        RevPtr<const RevVariable>                       rho;                                                                                    //!< The taxon sampling fraction
        RevPtr<const RevVariable>                       mu_rates;                                                                               //!< The extinction rates
        RevPtr<const RevVariable>                       mu_times;                                                                               //!< The extinction rate change times
        
    };
    
}

#endif
