#ifndef Dist_BDSTP_H
#define Dist_BDSTP_H

#include "RlBirthDeathProcess.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the Birth-Death-Sampling-Treatment Process
     *
     * The RevLanguage wrapper of the episodic or constant-rate birth-death-sampling-treatment process connects
     * the variables/parameters of the process and creates the internal BirthDeathSamplingTreatmentProcess object.
     * Please read the BirthDeathSamplingTreatmentProcess.h for more info.
     */
    class Dist_BDSTP : public BirthDeathProcess {

    public:
        Dist_BDSTP( void );

        // Basic utility functions
        Dist_BDSTP*                                             clone(void) const;                                                                      //!< Clone the object
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
        virtual void                                            addCommonRules(MemberRules&) const;                                                     //!< Add argument rules common to all BDSTP versions (default, FBD, phylodynamic)
        virtual void                                            addBurstRules(MemberRules &dist_member_rules) const;                                    //!< Add argument rules for burst/mass extinction events
        virtual RevBayesCore::DagNode*                          getRemovalProbability( void ) const;                                                    //!< Get removal probability

    private:

        RevPtr<const RevVariable>                               lambda;                                                                                 //!< The speciation rate(s)
        RevPtr<const RevVariable>                               mu;                                                                                     //!< The extinction rate(s)
        RevPtr<const RevVariable>                               phi;                                                                                    //!< The serial sampling rate(s)
        RevPtr<const RevVariable>                               r;                                                                                      //!< The treatment/removal probability(s)
        RevPtr<const RevVariable>                               Lambda;                                                                                 //!< The burst speciation rate(s)
        RevPtr<const RevVariable>                               Mu;                                                                                     //!< The mass extinction rate(s)
        RevPtr<const RevVariable>                               Phi;                                                                                    //!< The event sampling rate(s)
        RevPtr<const RevVariable>                               r_event;                                                                                //!< The event treatment/removal probability(s) - R, uses same timeline as Phi
        RevPtr<const RevVariable>                               timeline;                                                                               //!< The interval change times
        RevPtr<const RevVariable>                               lambda_timeline;                                                                        //!< The speciation rate change times
        RevPtr<const RevVariable>                               mu_timeline;                                                                            //!< The extinction rate change times
        RevPtr<const RevVariable>                               phi_timeline;                                                                           //!< The serial sampling rate change times
        RevPtr<const RevVariable>                               r_timeline;                                                                             //!< The change times for the treatment/removal probabilities
        RevPtr<const RevVariable>                               Lambda_timeline;                                                                        //!< The burst times
        RevPtr<const RevVariable>                               Mu_timeline;                                                                            //!< The mass extinction times
        RevPtr<const RevVariable>                               Phi_timeline;                                                                           //!< The event sampling times
        std::string                                             start_condition;                                                                        //!< The start condition of the process (rootAge/originAge)
        RevPtr<const RevVariable>                               initial_tree;                                                                           //!< User-specified initial tree (optional)
        RevPtr<const RevVariable>                               age_check_precision;                                                                    //!< Number of decimal places to use when checking the initial tree against taxon ages


    };

}

#endif
