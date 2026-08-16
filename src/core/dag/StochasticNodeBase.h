#ifndef StochasticNodeBase_H
#define StochasticNodeBase_H

#include "RbOrderedSet.h"
#include "SimulationConditions.h"

#include <cstddef>
#include <optional>
#include <vector>

namespace RevBayesCore {

    class DagNode;
    class Distribution;
    class DynamicNodeBase;

    // Put StochasticNode-specific stuff that doesn't depend on valueType here.
    class StochasticNodeBase {

    public:
        explicit StochasticNodeBase(Distribution *d);
        StochasticNodeBase(const StochasticNodeBase &n);
                                                           ~StochasticNodeBase(void);

        /**
         * Set directly the flag whether this node is clamped.
         * The caller needs to be responsible enough to know that we will assume
         * that the current value is the observed value.
         * We could use instead as well a call: clamp( getValue() );
         */
        void                                               setClamped(bool tf);

        /**
         * Unclamp this node. If I understand this correctly, we
         * can just set the clamped flag to false if we keep the
         * current value. We do not need to tell anyone that we
         * have changed value because we really haven't.
         */
        void                                               unclamp(void);

    protected:
        void                                               assign(const StochasticNodeBase &n, DagNode &owner);
        void                                               attachToDistributionParameters(DagNode &owner) const;
        void                                               bootstrap(DagNode &owner);
        virtual double                                     computeRecursiveIntegratedLnProbability(RbOrderedSet<DagNode*> &integratedParents, size_t index);
        void                                               detachFromDistributionParameters(DagNode &owner) const;
        void                                               getAffected(DagNode &owner, RbOrderedSet<DagNode*> &affected, const DagNode *affecter);
        void                                               getIntegratedParents(const DagNode &owner, RbOrderedSet<DagNode*> &integratedParents) const;
        double                                             getLnProbability(DagNode &owner);
        std::vector<double>                                getMixtureLikelihoods(const DagNode &owner, bool useLog) const;
        std::vector<double>                                getMixtureProbabilities(void) const;
        size_t                                             getNumberOfMixtureElements(void) const;
        std::vector<const DagNode*>                        getParents(void) const;
        double                                             getPrevLnProbability(void) const;
        bool                                               isClamped(void) const;
        bool                                               isIgnoredData(void) const;
        bool                                               isIntegratedOut(void) const;
        void                                               keepMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *affecter);
        void                                               redraw(DagNode &owner, SimulationCondition condition);
        void                                               reInitializeMe(void);
        void                                               restoreMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *restorer);
        void                                               setActivePIDSpecialized(size_t activePid, size_t numProcesses);
        void                                               setIgnoreData(const DagNode &owner, bool tf);
        void                                               setIgnoreRedraw(bool tf);
        void                                               setIntegratedOut(bool tf);
        void                                               setMcmcMode(bool tf);
        void                                               swapParent(DagNode &owner, const DagNode *oldParent, const DagNode *newParent);
        void                                               touchMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *toucher, bool touchAll);

        // protected members
        bool                                               clamped;
        bool                                               ignore_data;                                                                //!< The PDF for this node is set to 1, removing the effects of the node. Only for clamped nodes with no children.
        bool                                               ignore_redraw;
        mutable bool                                       integrated_out;
        std::optional<double>                              lnProb;                                                                     //!< Current log probability, or empty if not computed.
        std::optional<std::optional<double>>               stored_ln_prob;                                                             //!< Previous log probability if (a) there is a previous state and (b) the log probability for it is computed.
        Distribution                                      *distribution;
    };

}

#endif
