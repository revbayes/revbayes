#ifndef AbstractFossilizedBirthDeathRangeProcess_H
#define AbstractFossilizedBirthDeathRangeProcess_H

#include <algorithm>
#include <memory>

#include "RbConstants.h"
#include "MemberObject.h"
#include "RbVector.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"


namespace RevBayesCore {

    /**
     * One taxon's stratigraphic range: the four ages the model samples, ordered present to past,
     * and the bounds its occurrence record fixes them within. It owns the ordering, so the density
     * asks rather than restating the constraint.
     */
    struct RangeEntry {

        double  death = 0.0;                                //!< d_i, where the range ends: an extinction, or the marginalization limit
        double  last  = 0.0;                                //!< tau_K, the youngest augmented occurrence
        double  first = 0.0;                                //!< tau_1, the oldest augmented occurrence
        double  birth = 0.0;                                //!< b_i, where the range separates from its ancestor

        //!< The extremes are order statistics of the record, so every occurrence bounds them:
        //!< tau_1 is at least each reported minimum, tau_K at most each reported maximum.
        double  first_min = 0.0;                            //!< the oldest reported minimum
        double  first_max = RbConstants::Double::inf;       //!< the oldest reported maximum
        double  last_min  = 0.0;                            //!< the youngest reported minimum
        double  last_max  = RbConstants::Double::inf;       //!< the youngest reported maximum

        //!< Fewer than two reported occurrences, so tau_1 and tau_K are one age rather than two.
        bool    singleton = true;

        //!< Ages never decrease into the past, and each extreme sits within the bins that reported it.
        bool    isOrdered(double present) const
                {
                    if ( !( birth > first && first >= last && last >= death && death >= present ) ) return false;

                    // the extremes are order statistics of the record, so every reported bin bounds
                    // them: tau_1 clears every minimum, tau_K sits under every maximum
                    return first >= first_min && first <= first_max
                           && last >= last_min && last <= last_max;
                }

    };


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
     * of updateRanges(), which recomputes the origination and extinction times of each taxon
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
                                            const TypedDagNode<double>* origin = NULL,
                                            TypedDistribution<double>* origin_prior = NULL);  //!< Constructor

        virtual ~AbstractFossilizedBirthDeathRangeProcess(){};

        //!< The resampling move registers itself here, so the model can tell when it is missing
        void                                            setHasResampleMove(void) { has_resample_move = true; }

        //!< True when the augmented ages are elements of the owning distribution's value, so a
        //!< generic move on those elements samples them and restores them on rejection.
        virtual bool                                    augmentedAgesInValue(void) const { return false; }

        //!< A dnFossilRecord registers itself here. It carries the
        //!< whole sampling density, so without one psi has no data at all.
        void                                            setHasReportingNode(void) { has_reporting_node = true; }

        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;   //!< Expose the augmented first/last ages and the birth/death times for monitoring
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;              //!< Expose the origin, which no monitor can otherwise reach
        const std::vector<Taxon>&                       getTaxa() const { return taxa; }
        void                                            resampleFirstLast(size_t i);
        void                                            warnIfNoResampleMove(void) const;
        void                                            warnIfNoReportingNode(void) const;
        void                                            initializeFirstLast(size_t i);                               //!< Draw taxon i's augmented extremes nested (range_end <= last <= first) under the current reporting model.
        void                                            drawRanges();                                              //!< Draw an initial (range_start, range_end) and the augmented ages for every taxon. Shared by the matrix redraw and the tree (FBDSP) initial-value construction, which hangs a random budding topology on the ranges.
        double                                          computeLnFossilTotal();                                    //!< Total fossil-record log-density summed over taxa; used by a standalone dnFossilRecord node conditioned on this range process.
        //!< dnFossilRecord owns the reporting model and pushes it here, since the augmented ages
        //!< live on the range process. The cap is a constructor argument, so only a complete or
        //!< first/last record can be declared this way. Both augment the same ages, so the draw
        //!< the constructor already made stands.
        void                                            setCompleteRecord(bool comp) { record_complete = comp; }

        //!< Adopt a clamped record's occurrences as the data. The taxon set is the tree's tips, so
        //!< the caller must have checked that the names and order match; only the occurrences,
        //!< counts and extant status may differ.
        virtual void                                    setOccurrences(const std::vector<Taxon> &t);
        const std::vector<const DagNode*>&              getRangeParameters(void) const { return range_parameters; }   //!< The rate/timeline nodes a separate dnFossilRecord node adopts as parents.

    protected:
        virtual bool                                    marginalizesExtinction(void) const { return false; }    //!< True when range_end is integrated out, so a range ends at the marginalization limit and closes with p() rather than mu.
        virtual double                                  rangeEndTerm(size_t i, size_t di, double d) const { return log( death[di] ); }   //!< Log term closing range i at d > present.

        //!< Is tau_K an explicit latent? Not when unreported occurrences may lie below the youngest reported one.


        //!< Refresh the range table's state half from the value: the birth and death each process
        //!< keeps its ranges in, and for a tree the appearance ages the tips carry. The record
        //!< half (the occurrences and the bounds they fix) is updateRecord's, on the data's clock.
        //!< Re-establish what this process's own moves can break: the matrix's moves write one
        //!< column of a tied pair, and a tree's moves carry the ages its tips are authoritative for.
        virtual void                                    repairRanges(void) {}

        //!< Pull the ranges from the value. Runs on touch, where the proposal has already written
        //!< the value, so the density reads the table rather than refreshing it.
        virtual void                                    updateRanges() = 0;

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
        bool                                            record_complete;                                        //!< The declared reporting model: every occurrence reported, or the first/last rule
        double                                          max_present_age;                                        //!< Youngest occurrence maximum over all taxa; the timeline may not start above it

        void                                            updateRecord(void);                                //!< Refresh the range table's record half from the occurrences: the reported bins and the bounds they put on the appearance ages. Runs when the data changes, not per evaluation.

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

        std::vector<RangeEntry>                         ranges;                                                 //!< One per taxon: its four ages and the bounds its record fixes them within
        
        double                                          origin;                                                 //!< The origin time (oldest birth time)
        size_t                                          max_birth;                                              //!< Index of the taxon holding the oldest birth, refreshed by updateRanges

        // the following vectors are used internally for more efficient likelihood calculations and are filled by 'prepareProbComputation'
        mutable std::vector<double>                     birth;                                                  //!< The sorted speciation rates
        mutable std::vector<double>                     death;                                                  //!< The sorted extinction rates
        mutable std::vector<double>                     fossil;                                                 //!< The sorted fossil sampling rates
        mutable std::vector<double>                     times;                                                  //!< The sorted interval times
                        
        mutable std::vector<double>                     q_i;                                                    //!< Probability of no speciation or extinction event in each time interval
        mutable std::vector<double>                     q_tilde_i;                                              //!< Probability of no change in species identity in each time interval
        mutable std::vector<double>                     p_i;                                                    //!< Probability of leaving no sampled descendants from the end of each time interval
        mutable std::vector<double>                     pS_i;                                                   //!< Probability of leaving no descendants from the end of each time interval

                                
        //!< The one range mvResampleAugmentedAges drew, and the two ages it replaced. Taken in the
        //!< proposal, before the pull, so it undoes the proposal's own write and not the pull's.
        size_t                                          stored_range = 0;
        double                                          stored_first = 0.0;
        double                                          stored_last  = 0.0;

        //!< The table as the proposal found it. A rejected proposal restores the value, and the
        //!< table has to follow it back, or the next pull compares against the rejected state.
        std::vector<RangeEntry>                         stored_ranges;

        std::vector<double>                             partial_likelihood;                                     //!< Partial likelihood for each taxon
        std::vector<double>                             stored_likelihood;                                      //!< Stored partial likelihood for each taxon
                                
        std::vector<bool>                               dirty_taxa;                                             //!< Indicates whether partial likelihood needs updating
        
        bool                                            touched;                                               //!< Indicates whether any terms need updating
        bool                                            has_resample_move = false;              //!< Set by mvResampleAugmentedAges when it attaches
        bool                                            has_reporting_node = false;             //!< Set by dnFossilRecord when it attaches
        mutable bool                                    warned_no_resample = false;             //!< setMcmcMode fires more than once per run
        mutable bool                                    warned_no_reporting = false;
        bool                                            resampled;                                              //!< mvResampleAugmentedAges drew a new pair and the undo above is armed
    };
}

#endif
