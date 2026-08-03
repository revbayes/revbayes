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
     * @brief One taxon's stratigraphic range: the ages the model samples and the bounds its
     * occurrence record fixes them within.
     */
    struct RangeEntry {

        double  death = 0.0;                                //!< The extinction time, or the marginalization limit.
        double  last  = 0.0;                                //!< The youngest appearance, tau_K.
        double  first = 0.0;                                //!< The oldest appearance, tau_1.
        double  birth = 0.0;                                //!< The origination time.

        double  first_min = 0.0;                            //!< The oldest reported minimum.
        double  first_max = RbConstants::Double::inf;       //!< The oldest reported maximum.
        double  last_min  = 0.0;                            //!< The youngest reported minimum.
        double  last_max  = RbConstants::Double::inf;       //!< The youngest reported maximum.

        bool    singleton = true;                           //!< Fewer than two occurrences, so the appearances are one age.

        //!< Ages never decrease into the past, and each appearance sits within the bins reporting it.
        bool    isOrdered(double present) const
                {
                    if ( !( birth > first && first >= last && last >= death && death >= present ) ) return false;

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

        void                                            setHasResampleMove(void) { has_resample_move = true; }       //!< Registered by mvStratigraphicRange.
        void                                            setHasReportingNode(void) { has_reporting_node = true; }     //!< Registered by dnFossilRecord.

        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;   //!< Expose the appearances and the origination/extinction times.
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;              //!< Expose the origin.
        const std::vector<Taxon>&                       getTaxa() const { return taxa; }
        void                                            resampleFirstLast(size_t i);                                //!< Redraw taxon i's appearances within its reported bins.
        void                                            warnIfNoResampleMove(void) const;                           //!< Warn once when the appearances have no move.
        void                                            warnIfNoReportingNode(void) const;                          //!< Warn once when psi has no data.
        void                                            initializeFirstLast(size_t i);                              //!< Draw taxon i's appearances, nested within the range and its reported bins.
        void                                            drawRanges();                                               //!< Draw an initial range and appearances for every taxon.
        double                                          computeLnFossilTotal();                                     //!< Total fossil-record log-density, for a dnFossilRecord node.
        void                                            setCompleteRecord(bool comp) { record_complete = comp; }    //!< The reporting model, pushed here by dnFossilRecord.

        virtual void                                    setOccurrences(const std::vector<Taxon> &t);                //!< Adopt a clamped record's occurrences as the data.
        const std::vector<const DagNode*>&              getRangeParameters(void) const { return range_parameters; } //!< The rate and timeline nodes a dnFossilRecord adopts as parents.

    protected:
        virtual bool                                    marginalizesExtinction(void) const { return false; }        //!< True when the extinction times are integrated out.
        virtual double                                  rangeEndTerm(size_t i, size_t di, double d) const { return log( death[di] ); } //!< Log term closing range i at d.
        virtual double                                  rangeLnProb(size_t i);                                     //!< Birth-death density of range i. The one term a different birth-death model has to replace.
        virtual double                                  originLnProb(void);                                        //!< Terms attaching to the origin rather than to a range.
        virtual double                                  conditionLnProb(void) const;                               //!< The conditioning normalization, where it applies to the process as a whole.

        virtual void                                    repairRanges(void) {}                                       //!< Re-establish the invariants this process's own moves can break.
        virtual void                                    updateRanges() = 0;                                         //!< Pull the ranges from the value. Runs on touch, before the density reads them.

        virtual double                                  computeLnProbabilityRanges(bool force = false);
        double                                          computeLnFossilRecord(size_t i) const;                      //!< Fossil-record log-term for taxon i.

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
        bool                                            record_complete;                                        //!< Every occurrence reported, or the first/last rule.
        double                                          max_present_age;                                        //!< Youngest occurrence maximum; the timeline may not start above it.

        void                                            updateRecord(void);                                     //!< Refresh the reported bins and the bounds they put on the appearances.

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

        std::vector<RangeEntry>                         ranges;                                                 //!< One per taxon.
        
        double                                          origin;                                                 //!< The origin time (oldest birth time)
        size_t                                          max_birth;                                              //!< Index of the taxon holding the oldest origination time.

        // the following vectors are used internally for more efficient likelihood calculations and are filled by 'prepareProbComputation'
        mutable std::vector<double>                     birth;                                                  //!< The sorted speciation rates
        mutable std::vector<double>                     death;                                                  //!< The sorted extinction rates
        mutable std::vector<double>                     fossil;                                                 //!< The sorted fossil sampling rates
        mutable std::vector<double>                     times;                                                  //!< The sorted interval times
                        
        mutable std::vector<double>                     q_i;                                                    //!< Probability of no speciation or extinction event in each time interval
        mutable std::vector<double>                     q_tilde_i;                                              //!< Probability of no change in species identity in each time interval
        mutable std::vector<double>                     p_i;                                                    //!< Probability of leaving no sampled descendants from the end of each time interval
        mutable std::vector<double>                     pS_i;                                                   //!< Probability of leaving no descendants from the end of each time interval

                                
        size_t                                          stored_range = 0;                                       //!< The range mvStratigraphicRange drew, and the two ages it replaced.
        double                                          stored_first = 0.0;                                     //!< Taken in the proposal, before the pull, so it undoes the proposal's write.
        double                                          stored_last  = 0.0;

        std::vector<RangeEntry>                         stored_ranges;                                          //!< The table as the proposal found it, restored when it is rejected.

        std::vector<double>                             partial_likelihood;                                     //!< Partial likelihood for each taxon
        std::vector<double>                             stored_likelihood;                                      //!< Stored partial likelihood for each taxon
                                
        std::vector<bool>                               dirty_taxa;                                             //!< Indicates whether partial likelihood needs updating
        
        bool                                            touched;                                               //!< Indicates whether any terms need updating
        bool                                            has_resample_move = false;                              //!< Set by mvStratigraphicRange when it attaches.
        bool                                            has_reporting_node = false;                             //!< Set by dnFossilRecord when it attaches.
        mutable bool                                    warned_no_resample = false;                             //!< setMcmcMode fires more than once per run.
        mutable bool                                    warned_no_reporting = false;
        bool                                            resampled;                                              //!< The undo above is armed.
    };
}

#endif
