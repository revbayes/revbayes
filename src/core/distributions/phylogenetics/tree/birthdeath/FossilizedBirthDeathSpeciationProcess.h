#ifndef FossilizedBirthDeathSpeciationProcess_H
#define FossilizedBirthDeathSpeciationProcess_H

#include "RbVector.h"
#include "AbstractBirthDeathProcess.h"
#include "AbstractFossilizedBirthDeathRangeProcess.h"

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
     * @author The RevBayes Development Core Team (Walker Pett)
     * @since 2014-03-18, version 1.0
     *
     */
    class FossilizedBirthDeathSpeciationProcess : public AbstractBirthDeathProcess, public AbstractFossilizedBirthDeathRangeProcess {
        
        using AbstractBirthDeathProcess::taxa;

    public:
        FossilizedBirthDeathSpeciationProcess (const TypedDagNode<double>* root_age,
                                                      const DagNode *speciation,
                                                      const DagNode *extinction,
                                                      const DagNode *psi,
                                                      const DagNode *counts,
                                                      const TypedDagNode<double>* rho,
                                                      const TypedDagNode<RbVector<double> > *times,
                                                      const std::string &condition,
                                                      const std::vector<Taxon> &taxa,
                                                      bool use_origin,
                                                      bool presence_absence,
                                                      bool extended,
                                                      Tree *initial_tree,
                                                      long age_check_precision);  //!< Constructor
        
        // public member functions
        FossilizedBirthDeathSpeciationProcess*          clone(void) const override;                                               //!< Create an independent clone

        void                                            redrawValue(void) override;
        void                                            simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn) override;
        bool                                            allowsSA() override;                                                      //!< Checks if distribution is compatible with sampled ancestors

    protected:
        void                                            updateStartEndTimes() const override;
        int                                             updateStartEndTimes(const TopologyNode & ) const;

        double                                          pSurvival(double start, double end) const override;                       //!< Compute the probability of survival of the process (without incomplete taxon sampling).

        // Parameter management functions
        double                                          computeLnProbabilityTimes(void) const override;                           //!< Compute the log-transformed probability of the current value.
        double                                          computeLnProbabilityDivergenceTimes(void) const override;                 //!< Compute the log-transformed probability of the current value.

        double                                          lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const override { throw RbException("Cannot compute P(nTaxa)."); }
        double                                          lnProbTreeShape(void) const override;

        double                                          simulateDivergenceTime(double origin, double present) const override;     //!< Simulate a speciation event.
        std::vector<double>                             simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const override;                 //!< Simulate n speciation events.

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override; //!< Swap a parameter

    private:
        
        // helper functions
        double                                          getMaxTaxonAge( const TopologyNode& ) const;

        mutable std::vector<bool>                       I;
        bool                                            extended;
    };
}

#endif
