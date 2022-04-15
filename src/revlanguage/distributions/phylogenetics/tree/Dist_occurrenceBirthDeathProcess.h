#ifndef Dist_occurrenceBirthDeathProcess_H
#define Dist_occurrenceBirthDeathProcess_H

#include "RlBirthDeathProcess.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the piecewise constant Occurrence Birth Death Process
     *
     * The RevLanguage wrapper of the piecewise constant occurrence-birth-death process connects
     * the variables/parameters of the process and creates the internal OccurrenceBirthDeathProcess object.
     * Please read the OccurrenceBirthDeathProcess.h for more info.
     *
     *
     * @author Antoine Zwaans & Jérémy Andréoletti
     * @since 2020-03, version 1.0
     *
     */
    class Dist_occurrenceBirthDeathProcess : public BirthDeathProcess {

    public:
        Dist_occurrenceBirthDeathProcess( void );

        // Basic utility functions
        Dist_occurrenceBirthDeathProcess*                       clone(void) const;                              //!< Clone the object
        static const std::string&                               getClassType(void);                             //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                         //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;     //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;        //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                        //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                  //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::AbstractBirthDeathProcess*                createDistribution(void) const;

    protected:

        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);      //!< Set member variable


    private:

        RevPtr<const RevVariable>                               lambda;                                         //!< The speciation rate(s)
        RevPtr<const RevVariable>                               mu;                                             //!< The extinction rate(s)
        RevPtr<const RevVariable>                               psi;                                            //!< The serial sampling rate(s)
        RevPtr<const RevVariable>                               r;                                              //!< The probabilit(y|ies) of death upon sampling (treatment)
        RevPtr<const RevVariable>                               omega;                                          //!< The occurrence sampling rate(s)
        RevPtr<const RevVariable>                               rho;                                            //!< The present sampling rate
        RevPtr<const RevVariable>                               timeline;                                       //!< The interval change times
        RevPtr<const RevVariable>                               lambda_timeline;                                //!< The speciation/birth rate change times
        RevPtr<const RevVariable>                               mu_timeline;                                    //!< The extinction/death rate change times
        RevPtr<const RevVariable>                               psi_timeline;                                   //!< The serial sampling rate change times
        RevPtr<const RevVariable>                               r_timeline;                                     //!< The change times for the death upon sampling probabilities
        RevPtr<const RevVariable>                               omega_timeline;                                 //!< The occurrence sampling rate change times
        std::string                                             start_condition;                                //!< The start condition of the process (rootAge/originAge)
        RevPtr<const RevVariable>                               initial_tree;                                   //!< Optional initial tree
        RevPtr<const RevVariable>                               maxHiddenLin;                                   //!< The number of hidden lineages (algorithm accuracy)
        RevPtr<const RevVariable>                               occurrence_ages;                                //!< Occurrence ages
        RevPtr<const RevVariable>                               useMt;                                          //!< Forward traversal Mt algorithm (otherwise backward Lt)
        RevPtr<const RevVariable>                               verbose;                                        //!< Display warnings and information messages



    };

}

#endif
