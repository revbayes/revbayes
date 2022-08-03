#ifndef Dist_FBDP_H
#define Dist_FBDP_H

#include "RlBirthDeathProcess.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the constant-rate Fossilized-Birth-Death Process
     *
     * The RevLanguage wrapper of the episodic or constant-rate birth-death-sampling-treatment process connects
     * the variables/parameters of the process and creates the internal EpisodicBirthDeathSamplingTreatmentProcess object.
     * Please read the EpisodicBirthDeathSamplingTreatmentProcess.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *c
     */
    class Dist_FBDP : public BirthDeathProcess {

    public:
        Dist_FBDP( void );

        // Basic utility functions
        Dist_FBDP*                                              clone(void) const;                                                                      //!< Clone the object
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
        RevPtr<const RevVariable>                               phi;                                                                                    //!< The serial sampling rate(s)
        RevPtr<const RevVariable>                               r;                                                                                    //!< The taxon sampling fraction(s)
        RevPtr<const RevVariable>                               Lambda;                                                                                 //!< The speciation rate(s)
        RevPtr<const RevVariable>                               Mu;                                                                                     //!< The extinction rate(s)
        RevPtr<const RevVariable>                               Phi;                                                                                    //!< The serial sampling rate(s)
        RevPtr<const RevVariable>                               r_event;                                                                                    //!< The serial sampling rate(s)
        RevPtr<const RevVariable>                               timeline;                                                                               //!< The interval change times
        RevPtr<const RevVariable>                               lambda_timeline;                                                                        //!< The speciation rate change times
        RevPtr<const RevVariable>                               mu_timeline;                                                                            //!< The extinction rate change times
        RevPtr<const RevVariable>                               phi_timeline;                                                                           //!< The serial sampling rate change times
        RevPtr<const RevVariable>                               r_timeline;                                                                           //!< The change times for the death upon sampling probabilities
        RevPtr<const RevVariable>                               Lambda_timeline;                                                                        //!< The burst times
        RevPtr<const RevVariable>                               Mu_timeline;                                                                            //!< The mass extinction times
        RevPtr<const RevVariable>                               Phi_timeline;                                                                           //!< The event sampling times
        std::string                                             start_condition;                                                                     //!< The start condition of the process (rootAge/originAge)
        RevPtr<const RevVariable>                               initial_tree;

    };

}

#endif
