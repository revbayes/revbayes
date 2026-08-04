#ifndef BirthDeathWithRateshifts_H
#define BirthDeathWithRateshifts_H

#include "FossilizedBirthDeathRangeProcess.h"

namespace RevBayesCore {

    /**
     * @brief Birth-death-with-rateshifts (BDS) matrix distribution of Silvestro et al. (2019).
     *
     * The fossilized birth-death range process under complete lineage sampling: lineages are
     * treated as independent (no coexistence/gamma factor), and each range is normalized by the
     * fossil non-detection probability over its interval. Shares the range process's augmentation,
     * redraw and reporting machinery; only the likelihood differs.
     */
    class BirthDeathWithRateshifts : public FossilizedBirthDeathRangeProcess {

    public:
        BirthDeathWithRateshifts (const DagNode *speciation,
                                  const DagNode *extinction,
                                  const DagNode *psi,
                                  const TypedDagNode<double>* rho,
                                  const TypedDagNode<RbVector<double> > *times,
                                  const std::string &condition,
                                  const std::vector<Taxon> &taxa,
                                  bool complete_record,
                                  const TypedDagNode<double>* origin = NULL,
                                  double present = 0.0);

        BirthDeathWithRateshifts*                       clone(void) const override;

    protected:
        double                                          computeLnProbability(void) override;
        double                                          rangeLnProb(size_t i) override;                 //!< Waiting times and a Poisson record, with no p(t).
        double                                          originLnProb(void) override { return 0.0; }     //!< No origin term: this model has no p(t) to close the process with.
        double                                          conditionLnProb(void) const override { return 0.0; } //!< Conditioning on sampling is per range here, inside rangeLnProb.
    };
}

#endif
