#ifndef PiecewiseConstantFossilizedBirthDeathRangeProcess_H
#define PiecewiseConstantFossilizedBirthDeathRangeProcess_H

#include "AbstractPiecewiseConstantFossilizedRangeProcess.h"
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
     * @brief Piecewise-constant fossilized birth-death range distribution of origination extinction times matrix.
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
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class PiecewiseConstantFossilizedBirthDeathRangeProcess : public TypedDistribution<MatrixReal>, public AbstractPiecewiseConstantFossilizedRangeProcess {
        
    public:
        PiecewiseConstantFossilizedBirthDeathRangeProcess (const DagNode *speciation,
                                                      const DagNode *extinction,
                                                      const DagNode *psi,
                                                      const DagNode *counts,
                                                      const TypedDagNode<double>* rho,
                                                      const TypedDagNode<RbVector<double> > *times,
                                                      const std::string &condition,
                                                      const std::vector<Taxon> &taxa,
                                                      bool pa);  //!< Constructor
        
        // public member functions
        PiecewiseConstantFossilizedBirthDeathRangeProcess*   clone(void) const;                                         //!< Create an independent clone

    protected:
        void                                            updateStartEndTimes() const;

        // Parameter management functions
        double                                          computeLnProbability(void);                            //!< Compute the log-transformed probability of the current value.

        double                                          pSurvival(double start, double end) const;

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

        void                                            keepSpecialization(DagNode *toucher);
        void                                            restoreSpecialization(DagNode *toucher);
        void                                            touchSpecialization(DagNode *toucher, bool touchAll);

    private:
        
        // helper functions
        void                                            updateGamma(bool force = false);                             //!< Number of species alive at time t.
        void                                            redrawValue(void);

        mutable double                                  origin;

        std::vector<size_t>                             gamma_i;
        std::vector<std::vector<bool> >                 gamma_links;
        std::vector<bool>                               dirty_gamma;

        std::string                                     condition;
    };
}

#endif
