#ifndef Dist_OBDP2_H
#define Dist_OBDP2_H

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
    class Dist_OBDP2 : public BirthDeathProcess {

    public:
        Dist_OBDP2( void );

        // Basic utility functions
        Dist_OBDP2*                                             clone(void) const;                                                                      //!< Clone the object
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
        RevPtr<const RevVariable>                               r;
        RevPtr<const RevVariable>                               omega;                                                                                  //!< The occurrence sampling rate                                                                                 //!< The taxon sampling fraction(s)
        // RevPtr<const RevVariable>                               Lambda;                                                                                 //!< The speciation rate(s)
        // RevPtr<const RevVariable>                               Mu;                                                                                     //!< The extinction rate(s)
        RevPtr<const RevVariable>                               Phi;                                                                                    //!< The serial sampling rate(s)
        RevPtr<const RevVariable>                               timeline;                                                                               //!< The interval change times
        std::string                                             start_condition;                                                                        //!< The start condition of the process (rootAge/originAge)
        RevPtr<const RevVariable>                               initial_tree;                                                                           //!< Optional initial tree
        RevPtr<const RevVariable>                               maxHiddenLin;                                                                           //!< The number of hidden lineages (algorithm accuracy)
        RevPtr<const RevVariable>                               occurrence_ages;                                                                        //!< Occurrence ages
        RevPtr<const RevVariable>                               useMt;                                                                                  //!< Forward traversal Mt algorithm (otherwise backward Lt)
        RevPtr<const RevVariable>                               verbose;                                                                                //!< Display warnings and information messages



    };

}

#endif
