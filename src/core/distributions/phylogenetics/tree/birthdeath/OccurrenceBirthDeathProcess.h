#ifndef OccurrenceBirthDeathProcess_H
#define OccurrenceBirthDeathProcess_H

#include "AbstractBirthDeathProcess.h"
#include "RbVector.h"

#include <vector>
#include <set>

namespace RevBayesCore {

    class Clade;
    class Taxon;

    /**
    * @brief Piecewise-constant birth-death-sampling process with occurrences.
    *
    * The piecewise-constant occurrence birth-death process has constant rates for each time interval.
    * At the end of each interval there may be a rate-shift (jump for each of the rates).
    * Fossils with morphological characters are sampled with rate psi, fossil occurrences
    * without morphological at rate omega. Present lineages are conserved with probability rho.
    *
    * We assume that the rate vectors have one more element than the rate-change vectors.
    * Thus, one rate-change means always two interval, two rate-changes three interval, and so on.
    * The handling of heterogeneous/homogeneous rates
    * is adapted from the episodic birth-death-sampling-treatment process.
    *
    * @author Antoine Zwaans & Jérémy Andréoletti
    * @since 2020-03, version 1.0
    *
    */

    class OccurrenceBirthDeathProcess : public AbstractBirthDeathProcess {

    public:
        OccurrenceBirthDeathProcess (                           const TypedDagNode<double> *sa,
                                                                const DagNode *inspeciation,
                                                                const DagNode *inextinction,
                                                                const DagNode *inserialsampling,
                                                                const DagNode *intreatment,
                                                                const DagNode *inoccurrence,
                                                                const DagNode *ineventsampling,
                                                                const TypedDagNode< RbVector<double> > *speciation_timeline,
                                                                const TypedDagNode< RbVector<double> > *extinction_timeline,
                                                                const TypedDagNode< RbVector<double> > *sampling_timeline,
                                                                const TypedDagNode< RbVector<double> > *treatment_timeline,
                                                                const TypedDagNode< RbVector<double> > *occurence_timeline,
                                                                const TypedDagNode< RbVector<double> > *ht,
                                                                const std::string &cdt,
                                                                const std::vector<Taxon> &tn,
                                                                bool uo,
                                                                Tree *t,
                                                                const TypedDagNode<std::int64_t> *n,
                                                                const TypedDagNode< RevBayesCore::RbVector<double> > *O,
                                                                bool mt,
                                                                bool vb);                  //!< Constructor

        // public member functions
        OccurrenceBirthDeathProcess*   clone(void) const;        //!< Create an independent clone

    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);       //!< Swap a parameter

        // helper functions
        void                                                    addTimesToGlobalTimeline(std::set<double> &event_times, const TypedDagNode<RbVector<double> > *par_times) const;        //!< Adds timeline for parameter to set that we will use for global timeline
        void                                                    checkVectorSizes(const TypedDagNode<RbVector<double> >* v1, const TypedDagNode<RbVector<double> >* v2, int v1_minus_v2, const std::string& param_name, bool is_rate) const;
        void                                                    expandNonGlobalRateParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const; //!< Updates vector par such that it matches the global timeline
        void                                                    sortNonGlobalTimesAndVectorParameter(std::vector<double>& times, std::vector<double>& par) const;         //!< Sorts times to run from 0->inf, and orders par to match
        bool                                                    isConstantRate(void) const;                                             //!< Checks if we have a constant-rate process
        void                                                    sortGlobalTimesAndVectorParameter(void) const;                          //!< Sorts times to run from 0->inf, and orders ALL vector parameters to match

        double                                                  computeLnProbabilityDivergenceTimes(void) ;                            //!< Compute the log-transformed probability of the current value.
        double                                                  computeLnProbabilityTimes(void) const;                                 //!< Compute the log-transformed probability of the current value.
        size_t                                                  findIndex(double t) const;                                             //!< Find the index so that times[index-1] < t < times[index]
        size_t                                                  findIndex(double t, const std::vector<double>& timeline) const;

        void                                                    getOffset(void) const;
        double                                                  lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        void                                                    updateVectorParameters(void) const;
        double                                                  pSurvival(double start, double end) const;
        double                                                  simulateDivergenceTime(double origin, double present) const;           //!< Simulate a speciation event.
        void                                                    prepareTimeline(void) const;

        // members
        bool                                                    using_global_timeline;


        const TypedDagNode<double >*                            homogeneous_lambda;                         //!< The homogeneous birth rates.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_lambda;                       //!< The heterogeneous birth rates.
        const TypedDagNode<double >*                            homogeneous_mu;                             //!< The homogeneous death rates.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_mu;                           //!< The heterogeneous death rates.
        const TypedDagNode<double >*                            homogeneous_r;                              //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_r;                            //!< The heterogeneous conditional probability of death upon treatment.
        const TypedDagNode<double >*                            homogeneous_omega;                          //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_omega;                        //!< The heterogeneous conditional probability of death upon treatment.
        const TypedDagNode<double >*                            homogeneous_psi;                            //!< The homogeneous sampling rates.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_psi;                          //!< The heterogeneous sampling rates.
        const TypedDagNode<double >*                            homogeneous_rho;                            //!< The probability of sampling a tip at the present.

        const TypedDagNode<RbVector<double> >*                  interval_times_global;                      //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<RbVector<double> >*                  interval_times_speciation;                  //!< The user-specified non-zero times of rate shifts.
        const TypedDagNode<RbVector<double> >*                  interval_times_extinction;                  //!< The user-specified non-zero times of rate shifts.
        const TypedDagNode<RbVector<double> >*                  interval_times_sampling;                    //!< The user-specified non-zero times of rate shifts.
        const TypedDagNode<RbVector<double> >*                  interval_times_treatment;                   //!< The user-specified non-zero times of rate shifts.
        const TypedDagNode<RbVector<double> >*                  interval_times_occurrence;                  //!< The user-specified non-zero times of rate shifts.


        mutable std::vector<double>                             lambda_times;                               //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                             mu_times;                                   //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                             psi_times;                                  //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                             r_times;                                    //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                             omega_times;                                //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        mutable std::vector<double>                             global_timeline;                            //!< The times of the instantaneous events and rate shifts.

        const TypedDagNode<double >*                            start_age;                                  //!< The start age of the process
        const TypedDagNode<std::int64_t >*                              maxHiddenLin;                               //!< The maximal number of hidden lineages.
        const TypedDagNode< RevBayesCore::RbVector<double> >    *occurrence_ages;                           //!< Vector of occurrence times

        const bool                                              useMt;                                      //!< Forward traversal Mt algorithm (otherwise backward Lt)
        const std::string &                                     cond;                                       //!< Condition of the process (none/survival/survival2).
        const bool                                              verbose;                                    //!< Display warnings and information messages

        mutable std::vector<double>                             timeline;                                   //!< The times of the instantaneous events and rate shifts.
        mutable double                                          offset;                                     //!< In the case there the most recent tip is at time y, we internally adjust by this time and treat y as the present; this does not affect the boundary times of the rate shifts
        mutable std::vector<double>                             lambda;                                     //!< The speciation rates.
        mutable std::vector<double>                             mu;                                         //!< The extinction rates.
        mutable std::vector<double>                             psi;                                        //!< The fossil sampling rates.
        mutable std::vector<double>                             r;                                          //!< The removal probabilities after sampling.
        mutable std::vector<double>                             omega;                                      //!< The occurrence sampling rates.
        mutable std::vector<double>                             psi_event;                                  //!< The times of fossil sampling rate burst events.

        std::string spn = "speciation";
        std::string exn = "extinction";
        std::string smp = "sampling";
        std::string occ = "occurrence sampling";
        std::string trt = "(serial) treatment";
    };
}

#endif
