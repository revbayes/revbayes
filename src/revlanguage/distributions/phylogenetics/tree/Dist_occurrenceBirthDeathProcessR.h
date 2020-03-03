#ifndef Dist_OccurrenceBirthDeathProcessR_H
#define Dist_OccurrenceBirthDeathProcessR_H

#include "RlBirthDeathProcess.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the constant-rate Fossilized-Birth-Death Process
     *
     * The RevLanguage wrapper of the constant-rate fossilzed-birth-death process connects
     * the variables/parameters of the process and creates the internal ConstantRateFossilizedBirthDeathProcess object.
     * Please read the ConstantRateFossilizedBirthDeathProcess.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *c
     */
    class Dist_OccurrenceBirthDeathProcessR : public BirthDeathProcess {

    public:
        Dist_OccurrenceBirthDeathProcessR( void );

        // Basic utility functions
        Dist_OccurrenceBirthDeathProcessR*                      clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::AbstractBirthDeathProcess*                createDistribution(void) const;

    protected:

        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:

        RevPtr<const RevVariable>                               lambda;                                                                                 //!< The speciation rate(s)
        RevPtr<const RevVariable>                               mu;                                                                                     //!< The extinction rate(s)
        RevPtr<const RevVariable>                               psi;                                                                                    //!< The serial sampling rate(s)
        RevPtr<const RevVariable>                               rho;                                                                                    //!< The taxon sampling fraction(s)
        RevPtr<const RevVariable>                               timeline;                                                                               //!< The interval change times
        RevPtr<const RevVariable>                               lambda_timeline;                                                                        //!< The speciation rate change times
        RevPtr<const RevVariable>                               mu_timeline;                                                                            //!< The extinction rate change times
        RevPtr<const RevVariable>                               psi_timeline;                                                                           //!< The serial sampling rate change times
        RevPtr<const RevVariable>                               rho_timeline;                                                                           //!< The episodic taxon sampling fraction change times
        std::string                                             start_condition;                                                                     //!< The start condition of the process (rootAge/originAge)
        RevPtr<const RevVariable>                               initial_tree;

    };

}

#endif
