#ifndef PiecewiseConstantOccurrenceBirthDeathProcess_H
#define PiecewiseConstantOccurrenceBirthDeathProcess_H

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
    class PiecewiseConstantOccurrenceBirthDeathProcess : public AbstractBirthDeathProcess {

    public:
        PiecewiseConstantOccurrenceBirthDeathProcess (const TypedDagNode<double> *ra,
                                                      const TypedDagNode<double> *rho,
                                                      const DagNode *s,
                                                      const DagNode *e,
                                                      const DagNode *p,
                                                      const DagNode *o,
                                                      const DagNode *r,
                                                      const TypedDagNode< RbVector<double> > *ht,
                                                      const TypedDagNode< RbVector<double> > *lt,
                                                      const TypedDagNode< RbVector<double> > *mt,
                                                      const TypedDagNode< RbVector<double> > *pt,
                                                      const TypedDagNode< RbVector<double> > *ot,
                                                      const TypedDagNode< RbVector<double> > *rt,
                                                      const TypedDagNode<long> *n,
                                                      const std::string &cdt,
                                                      const std::vector<Taxon> &tn,
                                                      bool uo,
                                                      bool Mt,
                                                      TypedDagNode<Tree> *t);  //!< Constructor

        // public member functions
        PiecewiseConstantOccurrenceBirthDeathProcess*   clone(void) const;                                         //!< Create an independent clone

        double                                          getSpeciationRate( size_t index = 0 ) const;
        double                                          getExtinctionRate( size_t index = 0 ) const;
        double                                          getFossilSamplingRate( size_t index = 0 ) const;
        double                                          getOccurrenceSamplingRate( size_t index = 0 ) const;
        double                                          getRemovalProbability( size_t index = 0 ) const;

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
        const TypedDagNode< double > *                        start_age;                                       //!< Start age of the process.
        const TypedDagNode<double >*                    homogeneous_rho;                                       //!< The homogeneous extant sampling probs.
        const TypedDagNode<double >*                    homogeneous_lambda;                                    //!< The homogeneous speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                                  //!< The heterogeneous speciation rates.
        const TypedDagNode<double >*                    homogeneous_mu;                                        //!< The homogeneous extinction rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                                      //!< The heterogeneous extinction rates.
        const TypedDagNode<double >*                    homogeneous_psi;                                       //!< The homogeneous serial sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_psi;                                     //!< The heterogeneous serial sampling rates.
        const TypedDagNode<double >*                    homogeneous_omega;                                       //!< The homogeneous serial sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_omega;
        const TypedDagNode<double >*                    homogeneous_removalPr;                                       //!< The homogeneous serial sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_removalPr;


        const TypedDagNode<RbVector<double> >*          homogeneous_timeline;                                  //!< The times of the instantaneous rate change events.
        const TypedDagNode<RbVector<double> >*          lambda_timeline;                                       //!< The times of the instantaneous speciation rate change events.
        const TypedDagNode<RbVector<double> >*          mu_timeline;                                           //!< The times of the instantaneous extinction rate change events.
        const TypedDagNode<RbVector<double> >*          psi_timeline;                                          //!< The times of the instantaneous serial sampling rate change events.
        const TypedDagNode<RbVector<double> >*          omega_timeline;                                        //!< The times of the instantaneous sampling events.
        const TypedDagNode<RbVector<double> >*          removalPr_timeline;                                    //!< The times of the instantaneous sampling events.
        const TypedDagNode< long > *                    maxHiddenLin;                                          //!< The maximal number of hidden lineages.
        const bool                                      useMt;
        const std::string &                             cond;                                                  //!< Condition of the process ("time" or "survival")
        const bool                                      useOrigin;                                             //!< Start the process at the origin (otherwise root)

        mutable std::vector<double>                     timeline;
        mutable std::vector<double>                     lambda_times;
        mutable std::vector<double>                     mu_times;
        mutable std::vector<double>                     psi_times;
        mutable std::vector<double>                     omega_times;
        mutable std::vector<double>                     removalPr_times;
        mutable double                                  rho;
        mutable std::vector<double>                     lambda;
        mutable std::vector<double>                     mu;
        mutable std::vector<double>                     psi;
        mutable std::vector<double>                     omega;
        mutable std::vector<double>                     removalPr;


        bool                                            ascending;
    };
}

#endif
