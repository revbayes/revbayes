#ifndef Dist_GLHBDSP_H
#define Dist_GLHBDSP_H

#include "RlTimeTree.h"
#include "TreeDiscreteCharacterData.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the generalized lineage-heterogeneous birth-death-sampling process
     *
     * The RevLanguage wrapper of the generalized lineage-heterogeneous birth-death-sampling process connects
     * the variables/parameters of the process and creates the internal GeneralizedLineageHeterogeneousBirthDeathSamplingProcess object.
     * Please read the GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author Michael R. May and Xavier Meyer
     * @since 2020-03-09, version 1.0
     *c
     */

    class Dist_GLHBDSP : public TypedDistribution<TimeTree> {

    public:

    	Dist_GLHBDSP( void );
        virtual ~Dist_GLHBDSP(); //!< Virtual destructor

        // Basic utility functions
        Dist_GLHBDSP*                                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                     getDistributionMethods( void ) const;                                                   //!< Get the member methods
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::TypedDistribution<RevBayesCore::Tree>*    createDistribution(void) const;

    protected:

        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:

        // tolerances
        RevPtr<const RevVariable> abs_tol;          //!< The absolute tolerance of the integrator
        RevPtr<const RevVariable> rel_tol;          //!< The relative tolerance of the integrator

        // number of processors
        RevPtr<const RevVariable>  n_proc;          //!< The number of processors for parallel calculations

        // taxa
        RevPtr<const RevVariable>  taxa;            //!< The taxa in the tree

        // number of states
        RevPtr<const RevVariable>  n_states;        //!< The number of discrete states

        // age
        RevPtr<const RevVariable>  start_age;       //!< The age at the start of the process
        std::string                start_type;      //!< The start condition of the process (rootAge/originAge)

        // zero indexing
        RevPtr<const RevVariable>  zero_index;       //!< whether to zero index (default true)

        // regular events (and times)
        RevPtr<const RevVariable>  lambda;          //!< The speciation rates for each state for each epoch
        RevPtr<const RevVariable>  lambda_times;    //!< The times at which vectors of speciation-rates change
        RevPtr<const RevVariable>  mu;              //!< The extinction rates for each state for each epoch
        RevPtr<const RevVariable>  mu_times;        //!< The times at which vectors of extinction-rates change
        RevPtr<const RevVariable>  phi;             //!< The sampling rates for each state for each epoch
        RevPtr<const RevVariable>  phi_times;       //!< The times at which vectors of sampling-rates change
        RevPtr<const RevVariable>  delta;           //!< The destructive-sampling rates for each state for each epoch
        RevPtr<const RevVariable>  delta_times;     //!< The times at which vectors of destructive-sampling-rates change

        // mass events
        RevPtr<const RevVariable>  upsilon;         //!< The mass-speciation probabilities for each state at each time
        RevPtr<const RevVariable>  upsilon_times;   //!< The mass-speciation times
        RevPtr<const RevVariable>  gamma;           //!< The mass-extinction probabilities for each state at each time
        RevPtr<const RevVariable>  gamma_times;     //!< The mass-extinction times
        RevPtr<const RevVariable>  rho;             //!< The mass-sampling probabilities for each state at each time
        RevPtr<const RevVariable>  rho_times;       //!< The mass-sampling times
        RevPtr<const RevVariable>  xi;              //!< The mass-destructive-sampling probabilities for each state at each time
        RevPtr<const RevVariable>  xi_times;        //!< The mass-destructive-sampling times

        // state changes
        RevPtr<const RevVariable>  eta;             //!< The rates of change between each state for each epoch
        RevPtr<const RevVariable>  eta_times;       //!< The times at which the matrix of rates changes
        RevPtr<const RevVariable>  omega;           //!< The cladogenetic state-change probabilities for each epoch
        RevPtr<const RevVariable>  omega_times;     //!< The times at which the cladogenetic state-change probabilities change
        RevPtr<const RevVariable>  zeta;            //!< The probabilities of state change at each mass-extinction event

        // other parameters
        RevPtr<const RevVariable>  pi;               //!< The frequencies of the states at the root
        RevPtr<const RevVariable>  condition;        //!< The condition of the process

    };

}

#endif
