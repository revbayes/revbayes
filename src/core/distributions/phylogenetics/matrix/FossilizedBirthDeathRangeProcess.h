#ifndef FossilizedBirthDeathRangeProcess_H
#define FossilizedBirthDeathRangeProcess_H

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

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
    class FossilizedBirthDeathRangeProcess : public TypedDistribution<MatrixReal>, public AbstractFossilizedBirthDeathRangeProcess {
        
    public:
        FossilizedBirthDeathRangeProcess (const DagNode *speciation,
                                                      const DagNode *extinction,
                                                      const DagNode *psi,
                                                      const DagNode *counts,
                                                      const TypedDagNode<double>* rho,
                                                      const TypedDagNode<RbVector<double> > *times,
                                                      const std::string &condition,
                                                      const std::vector<Taxon> &taxa,
                                                      bool pa );  //!< Constructor
        
        // public member functions
        FossilizedBirthDeathRangeProcess*               clone(void) const override;                                                //!< Create an independent clone

    protected:
        void                                            updateStartEndTimes() const override;

        // Parameter management functions
        double                                          computeLnProbability(void) override;                                       //!< Compute the log-transformed probability of the current value.

        double                                          pSurvival(double start, double end) const;

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;  //!< Swap a parameter

        void                                            keepSpecialization(const DagNode *toucher) override;
        void                                            restoreSpecialization(const DagNode *toucher) override;
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll) override;

    private:
        
        // helper functions
        void                                            updateGamma(bool force = false);                                           //!< Number of species alive at time t.
        void                                            redrawValue(void) override;

        mutable double                                  origin;

        std::vector<size_t>                             gamma_i;
        std::vector<std::vector<bool> >                 gamma_links;
        std::vector<bool>                               dirty_gamma;

        std::string                                     condition;
    };
}

#endif
