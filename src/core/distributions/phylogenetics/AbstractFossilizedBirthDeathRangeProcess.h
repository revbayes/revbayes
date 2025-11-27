#ifndef AbstractFossilizedBirthDeathRangeProcess_H
#define AbstractFossilizedBirthDeathRangeProcess_H

#include "RbVector.h"
#include "Taxon.h"
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
    class AbstractFossilizedBirthDeathRangeProcess {
        
    public:
        AbstractFossilizedBirthDeathRangeProcess(const DagNode *speciation,
                                            const DagNode *extinction,
                                            const DagNode *psi,
                                            const TypedDagNode<double>* rho,
                                            const TypedDagNode<RbVector<double> > *times,
                                            const std::string &condition,
                                            const std::vector<Taxon> &taxa,
                                            bool complete,
                                            bool resampling);  //!< Constructor

        virtual ~AbstractFossilizedBirthDeathRangeProcess(){};

        std::vector<double>&                            getAges();
        void                                            resampleAge(size_t i);

    protected:
        virtual void                                    updateStartEndTimes() = 0;
        virtual double                                  computeLnProbabilityRanges(bool force = false);

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

        std::vector<const DagNode*>                     range_parameters;

        std::vector<double>                             b_i;                                                    //!< The birth times for each taxon
        std::vector<double>                             d_i;                                                    //!< The extinction times for each taxon
        std::vector<double>                             o_i;                                                    //!< The oldest minimum fossil age for each taxon
        std::vector<double>                             y_i;                                                    //!< The youngest maximum fossil age for each taxon
        
        double                                          origin;                                                 //!< The origin time (oldest birth time)

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
                                
        std::vector<double>                             age;                                                    //!< Age of the oldest occurence for each taxon
        std::vector<double>                             stored_age;                                             //!< Stored age of the oldest occurence for each taxon
                                
        std::vector<double>                             partial_likelihood;                                     //!< Partial likelihood for each taxon
        std::vector<double>                             stored_likelihood;                                      //!< Stored partial likelihood for each taxon
                                
        std::vector<bool>                               dirty_psi;                                              //!< Indicates whether fossil sampling terms need updating
        std::vector<bool>                               dirty_taxa;                                             //!< Indicates whether partial likelihood needs updating
        
        bool                                            complete;                                               //!< Indicates whether all fossil observations were included
        bool                                            touched;                                                //!< Indicates whether any terms need updating
        bool                                            resampled;                                              //!< Indicates whether any oldest occurrence ages were resampled
        bool                                            resampling;                                             //!< Indicates whether we are resampling oldest occurrence ages
    };
}

#endif
