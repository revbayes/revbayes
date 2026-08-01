#ifndef FossilizedBirthDeathSpeciationProcess_H
#define FossilizedBirthDeathSpeciationProcess_H

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "AbstractBirthDeathProcess.h"

namespace RevBayesCore {
    
    /**
     * @brief Piecewise-constant fossilized birth-death species distribution of extended trees.
     *
     * The piecewise-constant fossilized birth-death species process has constant rates for each time interval.
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
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2014-03-18, version 1.0
     *
     */
    class FossilizedBirthDeathSpeciationProcess : public AbstractBirthDeathProcess, public AbstractFossilizedBirthDeathRangeProcess {
        
        using AbstractFossilizedBirthDeathRangeProcess::taxa;

    public:
        FossilizedBirthDeathSpeciationProcess (const TypedDagNode<double>* ra,
                                      const DagNode *speciation,
                                      const DagNode *extinction,
                                      const DagNode *psi,
                                      const TypedDagNode<double>* rho,
                                      const DagNode *lambda_a,
                                      const DagNode *beta,
                                      const TypedDagNode<RbVector<double> > *times,
                                      const std::string &condition,
                                      const std::vector<Taxon> &taxa,
                                      bool complete_record,
                                      bool extended = true);  //!< Constructor
        
        // public member functions
        FossilizedBirthDeathSpeciationProcess*          clone(void) const override;                                //!< Create an independent clone

        void                                            setValue(Tree *v, bool force=false) override;                       //!< Clamping replaces the ranges the augmented ages were drawn against, so re-clip them
        void                                            redrawValue(void) override;
        void                                            redrawValue(SimulationCondition c) override;                        //!< The framework redraws through this overload, which must not reach the inherited simulator
        bool                                            redrawTopology(void);                                               //!< Redraw the budding topology uniformly with the ranges held fixed. Every compatible tree has one density, so this is a Gibbs step. Pure budding only, see hasAnagenesis() and hasSymmetricSpeciation().
        bool                                            hasAnagenesis(void) const;                                          //!< True when any interval's anagenetic rate is positive, which breaks the equal-density premise a budding topology Gibbs draw rests on.
        bool                                            hasSymmetricSpeciation(void) const;                                 //!< True when any interval's symmetric speciation probability is positive. Whether a given node may speciate symmetrically is symmetricAt().
        void                                            simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn) override;
        bool                                            allowsSampledAncestors(void) const override { return true; }                            //!< A sampled ancestor is an anagenetic speciation; with lambda_a = 0 its density is zero, so the flag stays on and the density does the rejecting.

    protected:
        bool                                            marginalizesExtinction(void) const override { return extended == false; }
        double                                          rangeEndTerm(size_t i, size_t di, double d) const override;

        void                                            updateRanges(void) override;
        void                                            setOccurrences(const std::vector<Taxon> &t) override;    //!< Sync the tree-side taxon copy, then adopt as the base does.
        double                                          symmetricAt(double age) const;                           //!< beta in the interval containing age, read from the parameter rather than the prepared cache, which the simulator paths run without.
        void                                            normalizeContinuationFlags(void);                        //!< Repair the whole tree. Refreshes the interval cache first, since the legal repair depends on beta at each node.
        void                                            normalizeContinuationFlags(const TopologyNode &node);    //!< Make each node name exactly one continuing child. Construction only: the density must reject an invalid state, not repair it.
        //!< What a subtree reports upward: the species running through it, or species == -1 when
        //!< that species ended below by symmetric speciation and the sampled ancestor above has yet
        //!< to name it. Invalidity travels separately, in invalid_continuation.
        struct RangeFlow { int species; double end_age; };
        RangeFlow                                       updateRanges(const TopologyNode & );

        double                                          pSurvival(double start, double end) const override;             //!< Compute the probability of survival of the process (without incomplete taxon sampling).

        // Parameter management functions
        double                                          computeLnProbabilityTimes(void) const override;                            //!< Compute the log-transformed probability of the current value.
        double                                          computeLnProbabilityDivergenceTimes(void) const override;            //!< Compute the log-transformed probability of the current value.

        //!< A non-extended tip is the augmented youngest age and must stay in its bin, but only for
        //!< an extinct taxon: an extant tip is pinned at the present, outside its fossil range.
        bool                                            tipAgeConstrainedToRange(const Taxon &t) const override { return extended == false && t.isExtinct(); }
        //!< The per-taxon check in the density covers this, and it alone knows which tips are exempt.
        bool                                            validatesTipAgesOnSet(void) const override { return false; }

        void                                            setMcmcMode(bool tf) override;

        bool                                            isExtended(void) const override { return extended; }            //!< An extended tree ends each range at the extinction time, so a tip may fall below its fossil age range. A non-extended tree ends it at the marginalization limit instead.

        double                                          lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const override { throw RbException("Cannot compute P(nTaxa)."); }
        double                                          lnProbTreeShape(void) const override;

        double                                          q(size_t i, double t, bool tilde = false) const override;

        double                                          simulateDivergenceTime(double origin, double present) const override;    //!< Simulate a speciation event.
        std::vector<double>                             simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const override;                 //!< Simulate n speciation events.

        void                                            keepSpecialization(const DagNode *toucher) override;
        void                                            restoreSpecialization(const DagNode *toucher) override;
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll) override;

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;                //!< Swap a parameter

        void                                            prepareProbComputation(void) const override;

    private:
        
        // helper functions
        double                                          getMaxTaxonAge( const TopologyNode& ) const;

        bool                                            extended;                                                //!< Tips are extinction times. When false the extinction times are marginalized out and each range ends at min(tau_K, youngest child birth).

        mutable std::vector<bool>                       I;                                                       //!< Indicates for each taxon whether the parent species was a sampled ancestor.
        mutable std::vector<bool>                       is_sa;                                                   //!< Indicates for each taxon whether its own range ends at a sampled ancestor node, so its lineage carries on.
        mutable std::vector<bool>                       ends_symmetric;                                          //!< Range i ends by symmetric speciation, so it is a speciation event and not an extinction.
        mutable double                                  budding_lnProb;                                          //!< Sum of log(1-beta) over budding events, accumulated by the tree pass.
        mutable bool                                    invalid_continuation = false;                            //!< Some node's children do not name exactly one continuation of its species.

        mutable std::vector<double>                     anagenetic;                                              //!< The sorted anagenetic speciation rates.
        mutable std::vector<double>                     symmetric;                                               //!< The sorted symmetric speciation probabilities.

        const TypedDagNode<double >*                    homogeneous_lambda_a;                                    //!< The homogeneous anagenetic speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda_a;                                  //!< The heterogeneous anagenetic speciation rates.
        const TypedDagNode<double >*                    homogeneous_beta;                                        //!< The homogeneous symmetric speciation probability.
        const TypedDagNode<RbVector<double> >*          heterogeneous_beta;                                      //!< The heterogeneous symmetric speciation probabilities.

    };
}

#endif
