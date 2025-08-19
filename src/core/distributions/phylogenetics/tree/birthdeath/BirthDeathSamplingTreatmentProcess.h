#ifndef BirthDeathSamplingTreatmentProcess_H
#define BirthDeathSamplingTreatmentProcess_H

#include "AbstractBirthDeathProcess.h"
#include "RbVector.h"

#include <vector>
#include <set>

namespace RevBayesCore {

    class Clade;
    class Taxon;

    /**
     * @brief Piecewise-constant birth-death-sampling-treatment process with mass extinctions, burst speciation, and event-sampling.
     *
     * The piecewise-constant birth-death process has constant rates for each time interval.
     * At the end of each time interval there may be an abrupt rate-shift (jump) for each
     * of the rates. Additionally, there may be sampling at the end of each interval.
     * Finally, fossils are sampled with rate psi, the others (fossils and extant taxa) are
     * sampled at sampling times (including the present).
     *
     * We assume that the rate vectors have one more element than the rate-change vectors.
     * Thus, one rate-change means always two interval, two rate-changes three interval, and so on.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class BirthDeathSamplingTreatmentProcess : public AbstractBirthDeathProcess {

    public:
        BirthDeathSamplingTreatmentProcess (const TypedDagNode<double> *ra,
                                                    const DagNode *speciation,
                                                    const DagNode *extinction,
                                                    const DagNode *sampling,
                                                    const DagNode *treatment,
                                                    const DagNode *event_speciation,
                                                    const DagNode *event_extinction,
                                                    const DagNode *event_sampling,
                                                    const DagNode *event_treatment,
                                                    const TypedDagNode< RbVector<double> > *timeline,
                                                    const TypedDagNode< RbVector<double> > *speciation_timeline,
                                                    const TypedDagNode< RbVector<double> > *extinction_timeline,
                                                    const TypedDagNode< RbVector<double> > *sampling_timeline,
                                                    const TypedDagNode< RbVector<double> > *treatment_timeline,
                                                    const TypedDagNode< RbVector<double> > *event_speciation_timeline,
                                                    const TypedDagNode< RbVector<double> > *event_extinction_timeline,
                                                    const TypedDagNode< RbVector<double> > *event_sampling_timeline,
                                                    const std::string &cdt,
                                                    const std::vector<Taxon> &tn,
                                                    bool uo,
                                                    Tree *t,
                                                    std::int64_t age_check_precision);  //!< Constructor

        // public member functions
        BirthDeathSamplingTreatmentProcess*             clone(void) const override;                                             //!< Create an independent clone
        void                                            redrawValue(SimulationCondition c = SimulationCondition::MCMC) override;//!< Draw a new random value from the distribution
        bool                                            allowsSA() override;                                                    //!< Checks if distribution is compatible with sampled ancestors
        virtual void                                    setValue(Tree *v, bool f=false) override;                               //!< Sets value to specified tree

    protected:
        // Parameter management functions
        double                                          computeLnProbabilityDivergenceTimes(void) const override;               //!< Compute the log-transformed probability of the current value.
        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;//!< Swap a parameter

        // helper functions
        void                                            addTimesToGlobalTimeline(std::set<double> &event_times, const TypedDagNode<RbVector<double> > *par_times) const;        //!< Adds timeline for parameter to set that we will use for global timeline
        void                                            checkVectorSizes(const TypedDagNode<RbVector<double> >* v1, const TypedDagNode<RbVector<double> >* v2, int v1_minus_v2, const std::string& param_name, bool is_rate) const;
        double                                          computeLnProbabilityTimes(void) const override;                         //!< Compute the log-transformed probability of the current value.
        bool                                            countAllNodes(void) const;                                              //!< Count bifurcating nodes, count all heterochronous nodes as either phi- or Phi-sampled and as either sampled ancestors or sampled extinct tips
        double                                          lnD(size_t i, double t) const;                                          //!< Branch-segment probability at time t with index i, using pre-computed vectors
        double                                          E(size_t i, double t, bool computeSurvival = false) const;                                       //!< Extinction probability at time t with index i, using pre-computed vectors
        void                                            expandNonGlobalProbabilityParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const; //!< Updates vector par such that it matches the global timeline
        void                                            expandNonGlobalRateParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const; //!< Updates vector par such that it matches the global timeline
        size_t                                          findIndex(double t) const;                                              //!< Find the index so that times[index-1] < t < times[index]
        size_t                                          findIndex(double t, const std::vector<double>& timeline) const;
        void                                            getOffset(void) const;
        bool                                            isConstantRate(void) const;                                             //!< Checks if we have a constant-rate process
        double                                          lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const override { throw RbException("Cannot compute P(nTaxa)."); }
        double                                          lnProbTreeShape(void) const override;
        void                                            prepareTimeline(void) const;
        void                                            prepareProbComputation(void) const override;
        double                                          pSampling(double t) const;
        double                                          pSurvival(double start, double end) const override;
        double                                          simulateDivergenceTime(double origin, double present) const override;   //!< Simulate a speciation event.
        void                                            sortGlobalTimesAndVectorParameter(void) const;                          //!< Sorts times to run from 0->inf, and orders ALL vector parameters to match
        void                                            sortNonGlobalTimesAndVectorParameter(std::vector<double>& times, std::vector<double>& par) const;     //!< Sorts times to run from 0->inf, and orders par to match
        int                                             survivors(double t) const;                                              //!< Number of species alive at time t.
        int                                             whichIntervalTime(double t) const;                                      //!< If a time corresponds to an interval/event time, returns that interval, otherwise returns -1

        // members
        bool                                            using_global_timeline;

        const TypedDagNode<double >*                    homogeneous_lambda;                                    //!< The homogeneous birth rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                                  //!< The heterogeneous birth rates.
        const TypedDagNode<double >*                    homogeneous_mu;                                        //!< The homogeneous death rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                                      //!< The heterogeneous death rates.
        const TypedDagNode<double >*                    homogeneous_r;                                         //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*          heterogeneous_r;                                       //!< The heterogeneous conditional probability of death upon treatment.
        const TypedDagNode<double >*                    homogeneous_phi;                                       //!< The homogeneous sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_phi;                                     //!< The heterogeneous sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_Lambda;                                  //!< The heterogeneous birth burst rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_Mu;                                      //!< The heterogeneous death burst (mass extinction) probabilities.
        const TypedDagNode<double >*                    homogeneous_Phi;                                       //!< The probability of sampling a tip at the present.
        const TypedDagNode<RbVector<double> >*          heterogeneous_Phi;                                     //!< The probability of sampling individuals at set time intervals.
        const TypedDagNode<RbVector<double> >*          heterogeneous_R;                                       //!< The heterogeneous conditional probability of death upon treatment.

        const TypedDagNode<RbVector<double> >*          interval_times_global;                                 //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*          interval_times_speciation;                             //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*          interval_times_extinction;                             //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*          interval_times_sampling;                               //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*          interval_times_treatment;                              //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*          interval_times_event_speciation;                       //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*          interval_times_event_extinction;                       //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*          interval_times_event_sampling;                         //!< The user-specified non-zero times of the instantaneous events and rate shifts.

        mutable std::vector<double>                     lambda_times;                                          //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                     mu_times;                                              //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                     phi_times;                                             //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                     r_times;                                               //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                     lambda_event_times;                                    //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                     mu_event_times;                                        //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                     phi_event_times;                                       //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                     global_timeline;                                       //!< The times of the instantaneous events and rate shifts.

        mutable std::vector<double>                     serial_tip_ages;                                       //!< The ages of all sampled dead lineages sampled by rate-sampling
        mutable std::vector<double>                     serial_sampled_ancestor_ages;                          //!< The ages of all sampled ancestors sampled by rate-sampling
        mutable std::vector<std::vector<double> >       event_tip_ages;                                        //!< The ages of all sampled dead lineages sampled at sampling events
        mutable std::vector<std::vector<double> >       event_sampled_ancestor_ages;                           //!< The ages of all sampled ancestors sampled at sampling events
        mutable std::vector<double>                     serial_bifurcation_times;                              //!< The ages of all bifurcation events in the tree NOT at a burst event
        mutable std::vector<std::vector<double> >       event_bifurcation_times;                               //!< The ages of all bifurcation events in the tree at burst events

        mutable double                                  offset;                                                //!< In the case there the most recent tip is at time y, we internally adjust by this time and treat y as the present; this does not affect the boundary times of the rate shifts
        int                                             num_extant_taxa;

        mutable std::vector<double>                     lambda;
        mutable std::vector<double>                     mu;
        mutable std::vector<double>                     phi;
        mutable std::vector<double>                     r;
        mutable std::vector<double>                     lambda_event;
        mutable std::vector<double>                     mu_event;
        mutable std::vector<double>                     phi_event;
        mutable std::vector<double>                     r_event;

        mutable std::vector<double>                     A_i;                                                   //!< Helper values
        mutable std::vector<double>                     B_i;                                                   //!< Helper values
        mutable std::vector<double>                     C_i;                                                   //!< Helper values
        mutable std::vector<double>                     E_previous;                                            //!< The probability that a lineage at time t has no sampled descendants.
        mutable std::vector<double>                     lnD_previous;                                          //!< The probability of an observed lineage

        mutable std::vector<double>                     A_survival_i;                                          //!< Helper values
        mutable std::vector<double>                     B_survival_i;                                          //!< Helper values
        mutable std::vector<double>                     C_survival_i;                                          //!< Helper values
        mutable std::vector<double>                     E_survival_previous;                                   //!< The probability that a lineage at time t has no sampled EXTANT descendants.

        std::string spn = "speciation";
        std::string exn = "extinction";
        std::string smp = "sampling";
        std::string trt = "(serial) treatment";
        std::string etrt = "(event) treatment";

        std::int64_t age_check_precision;
        std::vector<Taxon> taxa;
    };
}

#endif
