#ifndef PiecewiseConstantSerialSampledBirthDeathProcess_H
#define PiecewiseConstantSerialSampledBirthDeathProcess_H

#include "AbstractBirthDeathProcess.h"
#include "RbVector.h"

#include <vector>
#include <set>

namespace RevBayesCore {

    class Clade;
    class Taxon;

    /**
     * @brief Piecewise-constant fossilized birth-death process with serially sampled fossils.
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
    class PiecewiseConstantSerialSampledBirthDeathProcess : public AbstractBirthDeathProcess {

    public:
        PiecewiseConstantSerialSampledBirthDeathProcess (const TypedDagNode<double> *ra,
                                                      const DagNode *s,
                                                      const DagNode *e,
                                                      const DagNode *p,
                                                      const DagNode *r,
                                                      const TypedDagNode< RbVector<double> > *ht,
                                                      const TypedDagNode< RbVector<double> > *lt,
                                                      const TypedDagNode< RbVector<double> > *mt,
                                                      const TypedDagNode< RbVector<double> > *pt,
                                                      const TypedDagNode< RbVector<double> > *rt,
                                                      const std::string &cdt,
                                                      const std::vector<Taxon> &tn,
                                                      bool uo,
                                                      Tree *t);  //!< Constructor

        // public member functions
        PiecewiseConstantSerialSampledBirthDeathProcess*   clone(void) const;                                         //!< Create an independent clone

        double                                          getExtinctionRate( size_t index = 0 ) const;
        double                                          getSerialSamplingRate( size_t index = 0 ) const;
        double                                          getSpeciationRate( size_t index = 0 ) const;
        double                                          getTaxonSamplingProbability( size_t index = 0, size_t num_intervals = 1 ) const;

    protected:
        // Parameter management functions
        double                                          computeLnProbabilityDivergenceTimes(void) const;                            //!< Compute the log-transformed probability of the current value.
        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

        // helper functions
        double                                          computeLnProbabilityTimes(void) const;                     //!< Compute the log-transformed probability of the current value.
        size_t                                          l(double t) const;                                         //!< Find the index so that times[index-1] < t < times[index]
        double                                          lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        double                                          lnProbTreeShape(void) const;
        double                                          pSurvival(double start, double end) const;                 //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double                                          p(size_t i, double t) const;
        void                                            prepareProbComputation(void) const;
        double                                          q(size_t i, double t) const;
        double                                          simulateDivergenceTime(double origin, double present) const;    //!< Simulate a speciation event.
        int                                             survivors(double t) const;                                 //!< Number of species alive at time t.

        // members
        const TypedDagNode<double >*                    homogeneous_lambda;                                    //!< The homogeneous speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                                  //!< The heterogeneous speciation rates.
        const TypedDagNode<double >*                    homogeneous_mu;                                        //!< The homogeneous extinction rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                                      //!< The heterogeneous extinction rates.
        const TypedDagNode<double >*                    homogeneous_psi;                                       //!< The homogeneous serial sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_psi;                                     //!< The heterogeneous serial sampling rates.
        const TypedDagNode<double >*                    homogeneous_rho;                                       //!< The homogeneous extant sampling probs.
        const TypedDagNode<RbVector<double> >*          heterogeneous_rho;                                     //!< The heterogeneous extant sampling probs.

        const TypedDagNode<RbVector<double> >*          homogeneous_timeline;                                  //!< The times of the instantaneous rate change events.

        const TypedDagNode<RbVector<double> >*          lambda_timeline;                                       //!< The times of the instantaneous speciation rate change events.
        const TypedDagNode<RbVector<double> >*          mu_timeline;                                           //!< The times of the instantaneous extinction rate change events.
        const TypedDagNode<RbVector<double> >*          psi_timeline;                                          //!< The times of the instantaneous serial sampling rate change events.
        const TypedDagNode<RbVector<double> >*          rho_timeline;                                          //!< The times of the instantaneous sampling events.

        mutable std::vector<double>                     timeline;
        mutable std::vector<double>                     lambda_times;
        mutable std::vector<double>                     mu_times;
        mutable std::vector<double>                     psi_times;
        mutable std::vector<double>                     rho_times;
        mutable std::vector<double>                     lambda;
        mutable std::vector<double>                     mu;
        mutable std::vector<double>                     psi;
        mutable std::vector<double>                     rho;

        bool                                            ascending;
    };
}

#endif
