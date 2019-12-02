#ifndef Dist_BirthDeathBurstProcess_H
#define Dist_BirthDeathBurstProcess_H

#include "BirthDeathBurstProcess.h"
#include "RlBirthDeathProcess.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the episodic birth-death process.
     *
     * The RevLanguage wrapper of the episodic birth-death process connects
     * the variables/parameters of the process and creates the internal Dist_BirthDeathBurstProcess object.
     * Please read the Dist_BirthDeathBurstProcess.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2016-02-08, version 1.0
     *
     */
    class Dist_BirthDeathBurstProcess : public BirthDeathProcess {
        
    public:
        Dist_BirthDeathBurstProcess( void );
        
        // Basic utility functions
        Dist_BirthDeathBurstProcess*                        clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                            getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                 getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::BirthDeathBurstProcess*               createDistribution(void) const;
        
    protected:
        
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                           lambda;                                                                                 //!< The speciation rate
        RevPtr<const RevVariable>                           mu;                                                                                     //!< The extinction rate
        RevPtr<const RevVariable>                           burst_time;                                                                             //!< The time of the burst event
        RevPtr<const RevVariable>                           beta;                                                                                   //!< The probability of a speciation event at a burst event
        RevPtr<const RevVariable>                           rho;                                                                                    //!< The taxon sampling fraction

    };
    
}

#endif
