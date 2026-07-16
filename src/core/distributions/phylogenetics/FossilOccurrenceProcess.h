#ifndef FossilOccurrenceProcess_H
#define FossilOccurrenceProcess_H

#include "TypedDistribution.h"
#include "RbVector.h"
#include "Taxon.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * @brief Fossil sampling / reporting model, conditioned on a birth-death-range SKELETON.
     *
     * Observation half of the factored FBD-range model (path A of the skeleton/reporting
     * split). It contributes the per-taxon fossil-occurrence likelihood -- exactly the
     * `Psi[i]` block currently fused inside
     * AbstractFossilizedBirthDeathRangeProcess::computeLnProbabilityRanges() (lines ~307-521)
     * -- but as an explicit, clampable node:
     *
     *     occ ~ dnFossilOccurrences( skeleton )     // skeleton = dnFBDRP now; dnFBDSP later
     *     occ.clamp( taxa )                         // taxon-vector data type (prototype)
     *
     * DATA (this node's value): the taxon vector, whose per-taxon occurrence intervals and
     * counts are the observations.
     *
     * CONDITIONED-ON latent state -- true birth/death (b_i,d_i), the augmented extreme ages
     * (tau1 = first[i], tau_last = last[i]), the piecewise fossil-recovery rate psi and its
     * timeline -- is ALL read from the skeleton, never taken as this node's own parameters.
     * That structurally forbids a psi/timeline mismatch between the two halves (the coherence
     * point): there is exactly one psi and one timeline, owned by the skeleton.
     *
     * In path A the tau augmentation and the reporting-model choice (complete | firstlast |
     * uniform) live on the skeleton; this node delegates the term to
     * skeleton->computeLnFossilTotal() (added in increment 1). redrawValue() is a forward
     * Poisson-thin simulator => real posterior-predictive / isolated SBC on the reporting kernel.
     *
     * CONTRACT: the skeleton must provide traceable species lineages (budding taxonomy).
     * Satisfied by dnFBDRP always; by dnFBDSP only in the symmetric-speciation-prob -> 0 limit.
     *
     * OPEN (type-bridge, for the dev discussion): the skeleton is referenced as its process
     * object (AbstractFossilizedBirthDeathRangeProcess), not its Rev value, so one reporting
     * node can sit on either a MatrixReal-valued (dnFBDRP) or a TimeTree-valued (dnFBDSP)
     * skeleton. The prototype grabs it by downcasting the skeleton stochastic node's
     * distribution; the clean long-term form is a shared abstract Rev type.
     */
    class FossilOccurrenceProcess : public TypedDistribution< RbVector<Taxon> > {

    public:
        // ctor takes the skeleton's stochastic node (its distribution is downcast to
        // AbstractFossilizedBirthDeathRangeProcess) + the reported taxa (clamped data).
        FossilOccurrenceProcess(const DagNode *skeleton, const std::vector<Taxon> &taxa);
        virtual ~FossilOccurrenceProcess() {}

        FossilOccurrenceProcess*    clone(void) const override;
        double                      computeLnProbability(void) override;  // = skeleton->computeLnFossilTotal()
        void                        redrawValue(void) override;           // TODO: forward Poisson-thin simulator

    protected:
        void                        swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;

    private:
        const DagNode*                                   skeleton_node;   // owns b,d,tau,psi,timeline + reporting term
        const AbstractFossilizedBirthDeathRangeProcess*  skeleton;        // downcast view for computeLnFossilTotal()
        std::vector<Taxon>                               taxa;            // reported occurrences (data)
    };
}

#endif
