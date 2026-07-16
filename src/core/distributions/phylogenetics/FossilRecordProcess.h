#ifndef FossilRecordProcess_H
#define FossilRecordProcess_H

#include "TypedDistribution.h"
#include "RbVector.h"
#include "Taxon.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * @brief Observed fossil record, conditioned on a birth-death-range SKELETON.
     *
     * The observation half of the factored FBD-range model (Rev name: dnFossilRecord).
     * It contributes the per-taxon fossil-occurrence likelihood -- the term factored out of
     * AbstractFossilizedBirthDeathRangeProcess::computeLnProbabilityRanges() into
     * computeLnFossilTotal() -- as an explicit, clampable node:
     *
     *     rec ~ dnFossilRecord( skeleton )       // skeleton = dnFBDRP (matrix) or dnFBDSP (tree)
     *     rec.clamp( taxa )                      // taxon-vector data type (prototype)
     *
     * DATA (this node's value): the taxon vector, whose per-taxon occurrence intervals + counts
     * are the observations.
     *
     * CONDITIONED-ON latent state (b_i, d_i, the augmented extremes tau1/tau_last, the piecewise
     * fossil-recovery rate psi and its timeline) is ALL read from the skeleton via the shared
     * base, never taken as this node's own parameters -- one psi, one timeline, owned by the
     * skeleton (the coherence requirement). computeLnProbability() == skeleton->computeLnFossilTotal().
     *
     * REPORTING MODEL (complete | firstlast | uniform): a variant of ONE observation process,
     * differing only in the retention assumption -- so it is an OPTION, not a family of dists,
     * and it is applied per-taxon (effectiveReporting: uniform downgrades to complete below the
     * cap K). OPEN: in path A the reporting model currently lives on the skeleton (its `sampling`
     * member); the target is to move it here as a `reporting=` arg on dnFossilRecord (the skeleton
     * density q/q~ does not depend on it -- only the reporting term does). The dating timeline is a
     * further arg here (see redrawValue).
     *
     * redrawValue(): forward Poisson-thin simulator (real posterior-predictive / isolated SBC on the
     * reporting kernel). Encoding note: unbounded oldest-bin occurrences need max_age = Inf, which
     * requires the TaxonReader strtod fix and is only coherent under uniform reporting (see tasks).
     *
     * CONTRACT: the skeleton must provide traceable species lineages (budding taxonomy) -- dnFBDRP
     * always; dnFBDSP only in the symmetric-speciation-probability -> 0 limit.
     *
     * TYPE-BRIDGE (for the dev discussion): the skeleton is referenced as its process object
     * (AbstractFossilizedBirthDeathRangeProcess), not its Rev value, so one reporting node sits on
     * either a MatrixReal-valued (dnFBDRP) or TimeTree-valued (dnFBDSP) skeleton. The prototype grabs
     * it by downcasting the skeleton stochastic node's distribution; the clean long-term form is a
     * shared abstract Rev type.
     */
    class FossilRecordProcess : public TypedDistribution< RbVector<Taxon> > {

    public:
        FossilRecordProcess(const DagNode *skeleton, const std::vector<Taxon> &taxa);
        virtual ~FossilRecordProcess() {}

        FossilRecordProcess*        clone(void) const override;
        double                      computeLnProbability(void) override;  // == skeleton->computeLnFossilTotal()
        void                        redrawValue(void) override;           // TODO: forward Poisson-thin simulator

    protected:
        void                        swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;

        void                        resolveSkeleton(void);                // downcast skeleton_node's distribution to the range base

    private:
        const DagNode*                                   skeleton_node;   // owns b,d,tau,psi,timeline + reporting term
        AbstractFossilizedBirthDeathRangeProcess*        skeleton;        // downcast view for computeLnFossilTotal()
        std::vector<Taxon>                               taxa;            // reported occurrences (data)
    };
}

#endif
