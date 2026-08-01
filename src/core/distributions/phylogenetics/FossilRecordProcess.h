#ifndef FossilRecordProcess_H
#define FossilRecordProcess_H

#include "TypedDistribution.h"
#include "RbVector.h"
#include "Taxon.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * @brief Observed fossil record, conditioned on a birth-death range process (Rev: dnFossilRecord).
     *
     *     rec ~ dnFossilRecord( ranges, complete=false )   // ranges = dnFBDRP (matrix) or dnFBDSP (tree)
     *
     * The observation half of the factored FBD-range model: the per-taxon fossil-occurrence term
     * factored out of AbstractFossilizedBirthDeathRangeProcess::computeLnProbabilityRanges() into
     * computeLnFossilTotal(). The conditioned-on latent state (range_start, range_end, the appearances,
     * psi and its timeline) and the occurrences are all read from the range process, so one psi and
     * one timeline are owned there; computeLnProbability() == ranges->computeLnFossilTotal().
     *
     * The reporting model lives on the range process, set here via setCompleteRecord so both
     * nodes agree on the tau1 support and the reporting term.
     */
    class FossilRecordProcess : public TypedDistribution< RbVector<Taxon> > {

    public:
        FossilRecordProcess(const DagNode *ranges, bool complete);
        virtual ~FossilRecordProcess() {}

        FossilRecordProcess*        clone(void) const override;
        double                      computeLnProbability(void) override;  // == ranges->computeLnFossilTotal()
        void                        redrawValue(void) override;           // TODO: forward Poisson-thin simulator
        void                        setValue(RbVector<Taxon> *v, bool force=false) override;  // clamp: validate the record matches the range process

    protected:
        void                        swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;

        void                        resolveRanges(void);                // downcast ranges_node's distribution to the range base
        static AbstractFossilizedBirthDeathRangeProcess* rangesOf(const DagNode *n);  // shared downcast (MatrixReal or Tree node)

    private:
        const DagNode*                                   ranges_node;   // owns b,d,tau,psi,timeline + reporting term
        AbstractFossilizedBirthDeathRangeProcess*        ranges;        // downcast view for computeLnFossilTotal()
        bool                                             complete;      // every sampled occurrence reported; false is first/last. Pushed onto the range process
        std::vector<Taxon>                               taxa;          // reported occurrences (data)
    };
}

#endif
