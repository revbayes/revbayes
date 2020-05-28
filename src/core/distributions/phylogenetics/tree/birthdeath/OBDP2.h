#ifndef OBDP2_H
#define OBDP2_H

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
    class OBDP2 : public AbstractBirthDeathProcess {

    public:
        OBDP2 (                                       const TypedDagNode<double> *ra,
                                                      const DagNode *inspeciation,
                                                      const DagNode *inextinction,
                                                      const DagNode *inserialsampling,
                                                      const DagNode *intreatment,
                                                      const DagNode *inoccurrence,
                                                      const DagNode *ineventspeciation,
                                                      const DagNode *ineventextinction,
                                                      const DagNode *ineventsampling,
                                                      const TypedDagNode< RbVector<double> > *ht,
                                                      const std::string &cdt,
                                                      const std::vector<Taxon> &tn,
                                                      bool uo,
                                                      TypedDagNode<Tree> *t,
                                                      const TypedDagNode<long> *n,
                                                      const TypedDagNode< RevBayesCore::RbVector<double> > *O,
                                                      bool mt,
                                                      bool vb);  //!< Constructor

        // public member functions
        OBDP2*   clone(void) const;                                         //!< Create an independent clone

    protected:
        // Parameter management functions
        double                                          computeLnProbabilityDivergenceTimes(void) const;                            //!< Compute the log-transformed probability of the current value.
        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

        // helper functions
        double                                          computeLnProbabilityTimes(void) const;                     //!< Compute the log-transformed probability of the current value.
        void                                            countAllNodes(void) const;                                 //!< Count bifurcating nodes, count all heterochronous nodes as either phi- or Phi-sampled and as either sampled ancestors or sampled extinct tips
        // double                                          D(size_t i, double t) const;
        double                                          lnD(size_t i, double t) const;                                     //!< Branch-segment probability at time t with index i, using pre-computed vectors
        double                                          E(size_t i, double t) const;                                       //!< Extinction probability at time t with index i, using pre-computed vectors
        double                                          modifiedE(int i, double t, size_t condition_type) const;           //!< Extinction and extinction-like probabilities, recursively computed
        size_t                                          findIndex(double t) const;                                         //!< Find the index so that times[index-1] < t < times[index]
        void                                            getOffset(void) const;
        double                                          lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        double                                          lnProbTreeShape(void) const;
        void                                            updateVectorParameters(void) const;
        void                                            prepareProbComputation(void) const;
        double                                          pSurvival(double start, double end) const;
        double                                          simulateDivergenceTime(double origin, double present) const;    //!< Simulate a speciation event.
        int                                             survivors(double t) const;                                 //!< Number of species alive at time t.
        int                                             whichIntervalTime(double t) const;                                //!< If a time corresponds to an interval/event time, returns that interval, otherwise returns -1

        // members
        const TypedDagNode<double >*                    homogeneous_lambda;                                    //!< The homogeneous birth rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                                  //!< The heterogeneous birth rates.
        // const TypedDagNode<double >*                    homogeneous_Lambda;                                    //!< The homogeneous birth burst rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_Lambda;                                  //!< The heterogeneous birth burst rates.
        const TypedDagNode<double >*                    homogeneous_mu;                                        //!< The homogeneous death rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                                      //!< The heterogeneous death rates.
        const TypedDagNode<double >*                    homogeneous_r;                                         //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*          heterogeneous_r;
        const TypedDagNode<double >*                    homogeneous_o;                                         //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*          heterogeneous_o;                                         //!< The heterogeneous conditional probability of death upon treatment.
        // const TypedDagNode<double >*                    homogeneous_Mu;                                        //!< The homogeneous death burst (mass extinction) probabilities.
        const TypedDagNode<RbVector<double> >*          heterogeneous_Mu;                                      //!< The heterogeneous death burst (mass extinction) probabilities.
        const TypedDagNode<double >*                    homogeneous_phi;                                       //!< The homogeneous sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_phi;                                     //!< The heterogeneous sampling rates.
        const TypedDagNode<double >*                    homogeneous_Phi;                                       //!< The probability of sampling a tip at the present.
        const TypedDagNode<RbVector<double> >*          heterogeneous_Phi;                                     //!< The probability of sampling individuals at set time intervals.

        const TypedDagNode<RbVector<double> >*          interval_times;                                        //!< The user-specified non-zero times of the instantaneous events and rate shifts.
        const TypedDagNode<double >*                    start_age;
        const TypedDagNode<long >*                      maxHiddenLin;                                          //!< start age of the process
        const TypedDagNode< RevBayesCore::RbVector<double> > *occurrence_ages;                           //!< Vector of occurrence times
        const bool                                            useMt;                                 //!< Forward traversal Mt algorithm (otherwise backward Lt)
        const bool                                            useOrigin;
        const std::string &                                   cond;
        const bool                                            verbose;                               //!< Display warnings and information messages

        mutable std::vector<double>                     timeline;                                              //!< The times of the instantaneous events and rate shifts.

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
        mutable std::vector<double>                     omega;
        mutable std::vector<double>                     lambda_event;
        mutable std::vector<double>                     mu_event;
        mutable std::vector<double>                     phi_event;

        mutable std::vector<double>                     A_i;                                                   //!< Helper values
        mutable std::vector<double>                     B_i;                                                   //!< Helper values
        mutable std::vector<double>                     C_i;                                                   //!< Helper values
        mutable std::vector<double>                     E_previous;                                            //!< The probability that a lineage at time t has no sampled descendants.
        mutable std::vector<double>                     lnD_previous;                                          //!< The probability of an observed lineage
    };
}

#endif
