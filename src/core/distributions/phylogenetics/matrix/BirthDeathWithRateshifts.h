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
                                  size_t truncate_at,
                                  const TypedDagNode<double>* origin = NULL,
                                  bool report_internally = true);

        BirthDeathWithRateshifts*                       clone(void) const override;

    protected:
        double                                          computeLnProbability(void) override;
    };
}

#endif
