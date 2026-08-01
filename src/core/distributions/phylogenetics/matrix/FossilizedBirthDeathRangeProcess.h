#ifndef FossilizedBirthDeathRangeProcess_H
#define FossilizedBirthDeathRangeProcess_H

#include <set>
#include <utility>
#include "AbstractFossilizedBirthDeathRangeProcess.h"

#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
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
    class FossilizedBirthDeathRangeProcess : public TypedDistribution<MatrixReal>, public AbstractFossilizedBirthDeathRangeProcess {
        
    public:
        FossilizedBirthDeathRangeProcess (const DagNode *speciation,
										  const DagNode *extinction,
										  const DagNode *psi,
										  const TypedDagNode<double>* rho,
										  const TypedDagNode<RbVector<double> > *times,
										  const std::string &condition,
										  const std::vector<Taxon> &taxa,
										  bool complete_record,
                                          const TypedDagNode<double>* origin = NULL,
                                          TypedDistribution<double>* origin_prior = NULL);
        
        // public member functions
        void                                            setMcmcMode(bool tf) override;

        FossilizedBirthDeathRangeProcess*               clone(void) const override;                                         //!< Create an independent clone

        // Re-clip the appearances when the matrix is set externally (clamp);
        // MCMC moves mutate the value in place and never come through here.
        void                                            setValue(MatrixReal *v, bool force = false) override;

    protected:
        void                                            updateRanges() override;

        // Parameter management functions
        double                                          computeLnProbability(void) override;                                //!< Compute the log-transformed probability of the current value.

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;  //!< Swap a parameter

        std::vector<std::pair<size_t,double> >          stored_repairs;                                                     //!< Elements the repair wrote that no move stored, undone on restore.
        void                                            repairRanges(void) override;                                        //!< tau_1 and tau_K are one quantity below two occurrences; hold the columns equal. No move to follow: tau_K takes tau_1.
        void                                            repairRanges(const std::set<size_t> &touched);                      //!< Follow whichever of the two columns the move wrote.
        void                                            keepSpecialization(const DagNode *toucher) override;
        void                                            restoreSpecialization(const DagNode *toucher) override;
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll) override;

    private:
        
        // helper functions
        void                                            updateGamma(bool force = false);                                    //!< Number of species alive at time t.
        void                                            redrawValue(void) override;

        std::vector<size_t>                             gamma_i;            //!< The number of coexisting lineages at each taxon birth time
        std::vector<std::vector<bool> >                 gamma_links;        //!< A boolean matrix indicating which taxa are coexiting at each birth time
        std::vector<bool>                               dirty_gamma;        //!< Indicates whether gamma needs updating
    };
}

#endif
