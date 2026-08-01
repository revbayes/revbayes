#ifndef FossilizedBirthDeathSpeciationProcess_H
#define FossilizedBirthDeathSpeciationProcess_H

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "AbstractBirthDeathProcess.h"

namespace RevBayesCore {

    /**
     * @brief What the topology says about one taxon's range: how it begins, and how it ends.
     */
    struct SpeciesEntry {

        bool parent_is_sa   = false;    //!< The range this one budded from was sampled at that node.
        bool ends_at_sa     = false;    //!< The range ends at a sampled ancestor, so its lineage carries on.
        bool ends_symmetric = false;    //!< The range ends by symmetric speciation rather than extinction.

        bool operator==(const SpeciesEntry &e) const
                {
                    return parent_is_sa == e.parent_is_sa && ends_at_sa == e.ends_at_sa
                           && ends_symmetric == e.ends_symmetric;
                }

    };

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

        void                                            setValue(Tree *v, bool force=false) override;                       //!< Adopt a tree's annotations, or draw the state it does not carry.

        std::string                                     getHiddenStateString(void) const override;                          //!< The tree annotated with the range labels and appearances a Newick drops.
        void                                            setHiddenStateFromString(const std::string &s) override;
        void                                            redrawValue(void) override;
        void                                            redrawValue(SimulationCondition c) override;                        //!< The framework redraws through this overload, which must not reach the inherited simulator
        bool                                            redrawTopology(void);                                               //!< Redraw the budding topology uniformly with the ranges held fixed. Every compatible tree has one density, so this is a Gibbs step. Pure budding only, see hasAnagenesis() and hasSymmetricSpeciation().
        bool                                            hasAnagenesis(void) const;                                          //!< True when any interval's anagenetic rate is positive, which breaks the equal-density premise a budding topology Gibbs draw rests on.
        bool                                            hasSymmetricSpeciation(void) const;                                 //!< True when any interval's symmetric speciation probability is positive. Whether a given node may speciate symmetrically is symmetricAt().
        void                                            simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn) override;
        bool                                            allowsSampledAncestors(void) const override { return true; }        //!< With lambda_a = 0 the density is zero, so the density rejects rather than the flag.

    protected:
        bool                                            marginalizesExtinction(void) const override { return extended == false; }
        double                                          rangeEndTerm(size_t i, size_t di, double d) const override;

        void                                            repairRanges(void) override;                             //!< Re-take from the tip what a tree move carries there.
        void                                            updateRanges(void) override;
        void                                            setOccurrences(const std::vector<Taxon> &t) override;    //!< Sync the tree-side taxon copy, then adopt as the base does.
        double                                          symmetricAt(double age) const;                           //!< beta in the interval containing age, read from the parameter.
        void                                            normalizeContinuationFlags(void);                        //!< Repair the whole tree, refreshing the interval cache first.
        void                                            normalizeContinuationFlags(const TopologyNode &node);    //!< Make each node name exactly one continuing child.
        int                                             continuingSpecies(const TopologyNode &node) const;       //!< The taxon whose range this node belongs to, or -1 where none is named.

        //!< What a subtree reports upward: the range running through it, or -1 when that range
        //!< ended below and has yet to be named. Invalidity travels in invalid_continuation.
        struct RangeFlow { int species; double end_age; };
        RangeFlow                                       updateRanges(const TopologyNode & );

        double                                          pSurvival(double start, double end) const override;             //!< Compute the probability of survival of the process (without incomplete taxon sampling).

        // Parameter management functions
        double                                          computeLnProbabilityTimes(void) const override;                            //!< Compute the log-transformed probability of the current value.
        double                                          computeLnProbabilityDivergenceTimes(void) const override;            //!< Compute the log-transformed probability of the current value.

        bool                                            tipAgeConstrainedToRange(const Taxon &t) const override { return extended == false && t.isExtinct(); } //!< A non-extended extinct tip is its youngest appearance.
        bool                                            validatesTipAgesOnSet(void) const override { return false; }        //!< The density's per-taxon check knows which tips are exempt.

        void                                            setMcmcMode(bool tf) override;

        bool                                            isExtended(void) const override { return extended; }     //!< An extended tree ends each range at an extinction time.

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

        void                                            labelSpecies(TopologyNode &node, int s) const;           //!< Stamp the range each node belongs to, top down.
        bool                                            adoptSpeciesLabels(void);                               //!< Take the continuation flags from a tree's range labels.
        bool                                            adoptAppearances(void);                                 //!< Take the appearances from the tips' FAD and LAD.

        bool                                            extended;                                                //!< Tips are extinction times, else those are marginalized out.

        std::vector<SpeciesEntry>                       species;                                                 //!< One per taxon, as the last pull left it.
        std::vector<SpeciesEntry>                       next_species;                                            //!< Where the tree pass builds the replacement.
        std::vector<SpeciesEntry>                       stored_species;                                          //!< The table as the proposal found it.

        mutable double                                  budding_lnProb;                                          //!< Sum of log(1-beta) over budding events.
        mutable bool                                    invalid_continuation = false;                            //!< Some node's children do not name exactly one continuation.

        mutable std::vector<double>                     anagenetic;                                              //!< The sorted anagenetic speciation rates.
        mutable std::vector<double>                     symmetric;                                               //!< The sorted symmetric speciation probabilities.

        const TypedDagNode<double >*                    homogeneous_lambda_a;                                    //!< The homogeneous anagenetic speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda_a;                                  //!< The heterogeneous anagenetic speciation rates.
        const TypedDagNode<double >*                    homogeneous_beta;                                        //!< The homogeneous symmetric speciation probability.
        const TypedDagNode<RbVector<double> >*          heterogeneous_beta;                                      //!< The heterogeneous symmetric speciation probabilities.

    };
}

#endif
