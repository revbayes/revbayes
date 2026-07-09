#ifndef Dist_EpisodicStateDependentSpeciationExtinctionFossilizationProcess_H
#define Dist_EpisodicStateDependentSpeciationExtinctionFossilizationProcess_H

#include "RlTimeTree.h"
#include "TreeDiscreteCharacterData.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the episodic state-dependent birth-death-sampling process
     *
     * The RevLanguage wrapper of the episodic state-dependent birth-death-sampling process connects
     * the variables/parameters of the process and creates the internal EpisodicStateDependentSpeciationExtinctionFossilization object.
     * Please read the EpisodicStateDependentSpeciationExtinctionFossilization.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author Sebastian Höhna
     * @since 2026-07-07, version 1.4
     *c
     */

    class Dist_EpisodicStateDependentSpeciationExtinctionFossilizationProcess : public TypedDistribution<TimeTree> {

    public:

        Dist_EpisodicStateDependentSpeciationExtinctionFossilizationProcess( void );
        virtual ~Dist_EpisodicStateDependentSpeciationExtinctionFossilizationProcess(); //!< Virtual destructor

        // Basic utility functions
        Dist_EpisodicStateDependentSpeciationExtinctionFossilizationProcess*    clone(void) const;                                                                      //!< Clone the object
        static const std::string&                                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                                     getDistributionMethods( void ) const;                                                   //!< Get the member methods
        const TypeSpec&                                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                                      getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::TypedDistribution<RevBayesCore::Tree>*                    createDistribution(void) const;

    protected:

        void                                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:

        // tolerances
        RevPtr<const RevVariable>                                               abs_tol;          //!< The absolute tolerance of the integrator
        RevPtr<const RevVariable>                                               rel_tol;          //!< The relative tolerance of the integrator
        RevPtr<const RevVariable>                                               age_check_precision;//!< Number of decimal places to use when checking the initial tree against taxon ages

        // number of dense integrator attempts
        RevPtr<const RevVariable>                                               num_step;         //!< The maximum number of steps before giving up allowed for dense integrators

        // age
        RevPtr<const RevVariable>                                               start_age;       //!< The age at the start of the process
        std::string                                                             start_type;      //!< The start condition of the process (rootAge/originAge)

        // regular events (and times)
        RevPtr<const RevVariable>                                               lambda;          //!< The speciation rates for each state for each epoch
        RevPtr<const RevVariable>                                               lambda_times;    //!< The times at which vectors of speciation-rates change
        RevPtr<const RevVariable>                                               mu;              //!< The extinction rates for each state for each epoch
        RevPtr<const RevVariable>                                               mu_times;        //!< The times at which vectors of extinction-rates change
        RevPtr<const RevVariable>                                               phi;             //!< The sampling rates for each state for each epoch
        RevPtr<const RevVariable>                                               phi_times;       //!< The times at which vectors of sampling-rates change

        // mass events
        RevPtr<const RevVariable>                                               gamma;           //!< The mass-extinction probabilities for each state at each time
        RevPtr<const RevVariable>                                               gamma_times;     //!< The mass-extinction times

        // state changes
        RevPtr<const RevVariable>                                               eta;             //!< The rates of change between each state for each epoch
        RevPtr<const RevVariable>                                               eta_times;       //!< The times at which the matrix of rates changes
        RevPtr<const RevVariable>                                               omega;           //!< The cladogenetic state-change probabilities for each epoch
        RevPtr<const RevVariable>                                               omega_times;     //!< The times at which the cladogenetic state-change probabilities change

        // other parameters
        RevPtr<const RevVariable>                                               rho;              //!< The tip sampling probabilities
        RevPtr<const RevVariable>                                               pi;               //!< The frequencies of the states at the root
        RevPtr<const RevVariable>                                               condition;        //!< The condition of the process
        RevPtr<const RevVariable>                                               initial_tree;     //!< User-specified initial tree (optional)

    };

}

#endif
