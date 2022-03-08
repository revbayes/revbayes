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
     * Fossils with morphological caracters are sampled with rate psi, fossil occurrences
     * without morphological at rate omega. Present lineages are conserved with probability rho.
     *
     * We assume that the rate vectors have one more element than the rate-change vectors.
     * Thus, one rate-change means always two interval, two rate-changes three interval, and so on.
     * This implementation does not deal with rates that change at different times. If provided as vectors,
     * all vectors must be of equal length, to match the provided timeline.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Antoine Zwaans & Jérémy Andréoletti)
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
                                                                const TypedDagNode< RbVector<double> > *ht,
                                                                const std::string &cdt,
                                                                const std::vector<Taxon> &tn,
                                                                bool uo,
                                                                Tree *t,
                                                                const TypedDagNode<long> *n,
                                                                const TypedDagNode< RevBayesCore::RbVector<double> > *O,
                                                                bool mt,
                                                                bool vb);                  //!< Constructor

        // public member functions
        OccurrenceBirthDeathProcess*   clone(void) const;        //!< Create an independent clone

    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);       //!< Swap a parameter

        // helper functions
        double                                                  computeLnProbabilityDivergenceTimes(void) ;                            //!< Compute the log-transformed probability of the current value.
        double                                                  computeLnProbabilityTimes(void) const;                                       //!< Compute the log-transformed probability of the current value.
        size_t                                                  findIndex(double t) const;                                             //!< Find the index so that times[index-1] < t < times[index]
        void                                                    getOffset(void) const;
        double                                                  lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        void                                                    updateVectorParameters(void) const;
        double                                                  pSurvival(double start, double end) const;
        double                                                  simulateDivergenceTime(double origin, double present) const;           //!< Simulate a speciation event.

        // members
        const TypedDagNode<double >*                            homogeneous_lambda;                       //!< The homogeneous birth rates.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_lambda;                     //!< The heterogeneous birth rates.
        const TypedDagNode<double >*                            homogeneous_mu;                           //!< The homogeneous death rates.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_mu;                         //!< The heterogeneous death rates.
        const TypedDagNode<double >*                            homogeneous_r;                            //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_r;
        const TypedDagNode<double >*                            homogeneous_o;                            //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_o;                          //!< The heterogeneous conditional probability of death upon treatment.
        const TypedDagNode<double >*                            homogeneous_psi;                          //!< The homogeneous sampling rates.
        const TypedDagNode<RbVector<double> >*                  heterogeneous_psi;                        //!< The heterogeneous sampling rates.
        const TypedDagNode<double >*                            homogeneous_rho;                          //!< The probability of sampling a tip at the present.
        const TypedDagNode<RbVector<double> >*                  interval_times;                           //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<double >*                            start_age;                                //!< The start age of the process
        const TypedDagNode<long >*                              maxHiddenLin;                             //!< The maximal number of hidden lineages.
        const TypedDagNode< RevBayesCore::RbVector<double> >    *occurrence_ages;                         //!< Vector of occurrence times

        const bool                                              useMt;                                    //!< Forward traversal Mt algorithm (otherwise backward Lt)
        const bool                                              useOrigin;                                //!< If true the start age is the origin time otherwise the root age of the process.
        const std::string &                                     cond;                                     //!< Condition of the process (none/survival/survival2).
        const bool                                              verbose;                                  //!< Display warnings and information messages

        mutable std::vector<double>                             timeline;                                 //!< The times of the instantaneous events and rate shifts.
        mutable double                                          offset;                                   //!< In the case there the most recent tip is at time y, we internally adjust by this time and treat y as the present; this does not affect the boundary times of the rate shifts
        mutable std::vector<double>                             lambda;                                   //!< The speciation rates.
        mutable std::vector<double>                             mu;                                       //!< The extinction rates.
        mutable std::vector<double>                             psi;                                      //!< The fossil sampling rates.
        mutable std::vector<double>                             r;                                        //!< The removal probabilities after sampling.
        mutable std::vector<double>                             omega;                                    //!< The occurrence sampling rates.
        mutable std::vector<double>                             psi_event;                                //!< The times of fossil sampling rate burst events.
    };
}

#endif
