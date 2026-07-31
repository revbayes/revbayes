#ifndef AbstractFossilizedBirthDeathRangeProcess_H
#define AbstractFossilizedBirthDeathRangeProcess_H

#include <memory>

#include "MemberObject.h"
#include "RbVector.h"
#include "Taxon.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"


namespace RevBayesCore {

    /**
     * @brief Abstract piecewise-constant fossilized birth-death range process.
     *
     * The piecewise-constant fossilized birth-death range process has constant rates for each time interval.
     * At the end of each time interval there may be an abrupt rate-shift (jump) for each
     * of the rates. Additionally, there may be sampling at the end of each interval.
     * Finally, fossils are sampled with rate psi, the others (fossils and extant taxa) are
     * sampled at sampling times (including the present).
     *
     * We assume that the rate vectors have one more element than the rate-change vectors.
     * Thus, one rate-change means always two interval, two rate-changes three interval, and so on.
     * 
     * This is the base class for the fossilized birth-death range process, as well as
     * the fossilized birth-death range matrix process. Each process provides a different implementation
     * of updateStartEndTimes(), which recomputes the origination and extinction times of each taxon
     * based on either a matrix or tree data structure.
     * Likelihood computations are otherwise the same for all range-based fossilized birth-death processes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class AbstractFossilizedBirthDeathRangeProcess : public MemberObject< RbVector<double> >, public MemberObject<double> {
        
    public:
        AbstractFossilizedBirthDeathRangeProcess(const DagNode *speciation,
                                            const DagNode *extinction,
                                            const DagNode *psi,
                                            const TypedDagNode<double>* rho,
                                            const TypedDagNode<RbVector<double> > *times,
                                            const std::string &condition,
                                            const std::vector<Taxon> &taxa,
                                            bool complete_record,
                                            size_t truncate_at,
                                            const TypedDagNode<double>* origin = NULL,
                                            TypedDistribution<double>* origin_prior = NULL);  //!< Constructor

        virtual ~AbstractFossilizedBirthDeathRangeProcess(){};

        //!< The resampling move registers itself here, so the model can tell when it is missing
        void                                            setHasResampleMove(void) { has_resample_move = true; }

        //!< True when the augmented ages are elements of the owning distribution's value, so a
        //!< generic move on those elements samples them and restores them on rejection.
        virtual bool                                    augmentedAgesInValue(void) const { return false; }

        //!< A dnFossilRecord registers itself here. With report_internally false it carries the
        //!< whole sampling density, so without one psi has no data at all.
        void                                            setHasReportingNode(void) { has_reporting_node = true; }

        std::vector<double>&                            getAges();
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;   //!< Expose the augmented first/last ages and the birth/death times for monitoring
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;              //!< Expose the origin, which no monitor can otherwise reach
        const std::vector<Taxon>&                       getTaxa() const { return taxa; }
        double                                          getPresent(void) const { return times.front(); }                                   //!< Age of the present
        double                                          getOrigin(void) const { return origin; }                                            //!< The origin of the process, which is the oldest birth
        void                                            resampleFirstLast(size_t i);
        void                                            warnIfNoResampleMove(void) const;
        void                                            warnIfNoReportingNode(void) const;
        void                                            drawAugmentedAges(size_t i);                               //!< Draw taxon i's augmented extremes nested (range_end <= last <= first) under the current reporting model.
        void                                            drawRanges();                                              //!< Draw an initial (range_start, range_end) and the augmented ages for every taxon. Shared by the matrix redraw and the tree (FBDSP) initial-value construction, which hangs a random budding topology on the ranges.
        double                                          computeLnFossilTotal();                                    //!< Total fossil-record log-density summed over taxa; used by a standalone dnFossilRecord node conditioned on this range process (self-contained: refreshes rate cache + start/end times).
        //!< dnFossilRecord owns the reporting model and pushes it here, since the augmented ages
        //!< live on the range process. The cap is a constructor argument, so only a complete or
        //!< first/last record can be declared this way. Both augment the same ages, so the draw
        //!< the constructor already made stands.
        void                                            setCompleteRecord(bool comp) { record_complete = comp; complete = std::vector<bool>(taxa.size(), comp); }

        //!< Adopt a clamped record's occurrences as the data. The taxon set is the tree's tips, so
        //!< the caller must have checked that the names and order match; only the occurrences,
        //!< counts and extant status may differ.
        void                                            setOccurrences(const std::vector<Taxon> &t);
        const std::vector<const DagNode*>&              getRangeParameters(void) const { return range_parameters; }   //!< The rate/timeline nodes a separate dnFossilRecord node adopts as parents.

    protected:
        virtual bool                                    marginalizesExtinction(void) const { return false; }    //!< True when range_end is integrated out, so a range ends at the marginalization limit and closes with p() rather than mu.
        virtual double                                  rangeEndTerm(size_t i, size_t di, double d) const { return log( death[di] ); }   //!< Log term closing range i at d > present.

        //!< Is tau_K an explicit latent? Not when unreported occurrences may lie below the youngest reported one.
        bool                                            augmentsYoungest(size_t i) const { return occurrence_counts[i] >= 2 && truncated[i] == false; }

        //!< Adapters to the owning distribution, which this mixin is not a base of: its density, and
        //!< a fresh draw of its value.
        virtual double                                  ownLnProbability() = 0;
        virtual void                                    ownRedrawValue() = 0;

        //!< The augmented extremes are latent state that no node's value holds, so a checkpoint
        //!< drops them. These move them into the owning distribution's value and back. No-ops until
        //!< each process has somewhere to put them (a wider matrix, a tree that carries them).
        virtual void                                    pushAugmentedToValue() {}
        virtual void                                    pullAugmentedFromValue() {}

        void                                            clipAugmentedAges();                                    //!< Put every augmented extreme back inside its bin.
        bool                                            startsFinite();                                         //!< Does the current value score finitely under both the range and the reporting term?
        //!< The value was drawn against the old occurrences, so move it back inside the new bins,
        //!< redrawing it if no clipping suffices. False when nothing found puts it in the support.
        virtual bool                                    reclipToOccurrences();

        virtual void                                    updateStartEndTimes() = 0;
        virtual double                                  computeLnProbabilityRanges(bool force = false);
        double                                          computeLnFossilRecord(size_t i) const;              //!< Fossil-record (occurrence) log-term for taxon i, factored out of computeLnProbabilityRanges (range/reporting split).

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

        // helper functions
        size_t                                          findIndex(double t) const;                             //!< Find the index so that times[index-1] < t < times[index]
        double                                          p(size_t i, double t, bool survival = false) const;
        virtual double                                  q(size_t i, double t, bool tilde = false) const;

        virtual void                                    prepareProbComputation(void) const;

        void                                            keepSpecialization(const DagNode *toucher);                  /* NOT VIRTUAL */
        void                                            restoreSpecialization(const DagNode *toucher);               /* NOT VIRTUAL */
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll);  /* NOT VIRTUAL */

        std::vector<Taxon>                              taxa;                                                  //!< Taxa that will be attached to new simulated trees.
        std::string                                     condition;
        size_t                                          record_cap;                                            //!< The reporting cap K, kept so the per-taxon record model can be re-derived on a clamp
        bool                                            record_complete;                                        //!< The declared reporting model, same
        double                                          max_present_age;                                        //!< Youngest occurrence maximum over all taxa; the timeline may not start above it

        void                                            deriveRecordModel(void);                                //!< Recompute counts, first_min/last_max and the per-taxon reporting model from the occurrences

        size_t                                          num_intervals;

        // members
        const TypedDagNode<double >*                    homogeneous_lambda;                                    //!< The homogeneous speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                                  //!< The heterogeneous speciation rates.
        const TypedDagNode<double >*                    homogeneous_mu;                                        //!< The homogeneous speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                                      //!< The heterogeneous speciation rates.
        const TypedDagNode<double >*                    homogeneous_psi;                                       //!< The homogeneous speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_psi;                                     //!< The heterogeneous speciation rates.
        const TypedDagNode<double >*                    homogeneous_rho;                                       //!< The homogeneous speciation rates.
        const TypedDagNode<RbVector<double> >*          timeline;                                              //!< The times of the instantaneous sampling events.
        const TypedDagNode<double >*                    origin_age;                                            //!< The optional origin time of the process (NULL = oldest sampled birth).
        std::shared_ptr<TypedDistribution<double> >     origin_prior;                                          //!< Optional prior evaluated at the oldest birth, which is the origin.

        std::vector<const DagNode*>                     range_parameters;

        std::vector<double>                             range_start;                                                    //!< The birth times for each taxon
        std::vector<double>                             range_end;                                                    //!< The extinction times for each taxon
        std::vector<double>                             first_min;                                                    //!< The oldest minimum fossil age for each taxon
        std::vector<double>                             last_max;                                                    //!< The youngest maximum fossil age for each taxon
        
        double                                          origin;                                                 //!< The origin time (oldest birth time)
        size_t                                          max_birth;                                              //!< Index of the taxon holding the oldest birth, refreshed by updateStartEndTimes

        // the following vectors are used internally for more efficient likelihood calculations and are filled by 'prepareProbComputation'
        mutable std::vector<double>                     birth;                                                  //!< The sorted speciation rates
        mutable std::vector<double>                     death;                                                  //!< The sorted extinction rates
        mutable std::vector<double>                     fossil;                                                 //!< The sorted fossil sampling rates
        mutable std::vector<double>                     times;                                                  //!< The sorted interval times
                        
        mutable std::vector<double>                     q_i;                                                    //!< Probability of no speciation or extinction event in each time interval
        mutable std::vector<double>                     q_tilde_i;                                              //!< Probability of no change in species identity in each time interval
        mutable std::vector<double>                     p_i;                                                    //!< Probability of leaving no sampled descendants from the end of each time interval
        mutable std::vector<double>                     pS_i;                                                   //!< Probability of leaving no descendants from the end of each time interval

        std::vector<double>                             Psi;                                                    //!< Fossil sampling terms computed for each taxon
        std::vector<double>                             stored_Psi;                                             //!< Stored fossil sampling terms
                                
        std::vector<double>                             first;                                                  //!< Age of the oldest occurence for each taxon
        std::vector<double>                             stored_first;                                           //!< Stored age of the oldest occurence for each taxon
        std::vector<double>                             last;                                                   //!< Augmented age of the youngest occurrence for each taxon (first/last conditioning)
        std::vector<double>                             stored_last;                                            //!< Stored augmented age of the youngest occurrence for each taxon
                                
        std::vector<double>                             partial_likelihood;                                     //!< Partial likelihood for each taxon
        std::vector<double>                             stored_likelihood;                                      //!< Stored partial likelihood for each taxon
                                
        std::vector<bool>                               dirty_psi;                                              //!< Indicates whether fossil sampling terms need updating
        std::vector<bool>                               dirty_taxa;                                             //!< Indicates whether partial likelihood needs updating
        
        std::vector<size_t>                             occurrence_counts;                                      //!< Number of reported occurrences for each taxon
        std::vector<bool>                               truncated;                                              //!< Taxa reported up to the cap K, whose unreported specimens are marginalized
        std::vector<bool>                               complete;                                               //!< Taxa whose whole record is reported. Neither flag set is the first/last rule: both extremes reported, the interior count marginalized
        bool                                            touched;                                               //!< Indicates whether any terms need updating
        bool                                            has_resample_move = false;              //!< Set by mvResampleAugmentedAges when it attaches
        bool                                            has_reporting_node = false;             //!< Set by dnFossilRecord when it attaches
        mutable bool                                    warned_no_resample = false;             //!< setMcmcMode fires more than once per run
        mutable bool                                    warned_no_reporting = false;
        bool                                            resampled;                                              //!< Indicates whether any oldest occurrence ages were resampled
        bool                                            report_internally;                                      //!< If true (default) computeLnProbabilityRanges adds the reporting term inline (fused facade); if false the term is omitted (bare range process) and supplied by a separate dnFossilRecord node.
    };
}

#endif
