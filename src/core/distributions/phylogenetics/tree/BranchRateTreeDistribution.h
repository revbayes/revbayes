#ifndef BranchRateTreeDistribution_H
#define BranchRateTreeDistribution_H

#include <stddef.h>
#include <vector>

#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"
#include "Clade.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;


    class BranchRateTreeDistribution : public TypedDistribution<Tree>, public TreeChangeEventListener {

    public:
        BranchRateTreeDistribution(const TypedDagNode<Tree>* tt, TypedDistribution<double>* brp, TypedDagNode<double> *rbf);
        BranchRateTreeDistribution(const BranchRateTreeDistribution &d);
        virtual                                            ~BranchRateTreeDistribution(void);                                                                    //!< Virtual destructor

        BranchRateTreeDistribution&                         operator=(const BranchRateTreeDistribution &d);

        // public member functions
        BranchRateTreeDistribution*                         clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        virtual void                                        fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                 //!< This node was changed in the tree
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                        //!< Set the current value, e.g. attach an observation (clamp)

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter

    private:

        // helper functions
        RbBitSet                                            collectTreeSample(const TopologyNode& n, RbBitSet& in, std::map<RbBitSet, double>& bl);
        RbBitSet                                            collectSplits(const TopologyNode& n, RbBitSet& in, std::vector<RbBitSet>& s) const;
        void                                                simulateTree(void);
        void                                                simulateClade(std::vector<TopologyNode*> &n);                           //!< Simulate n speciation events.

        // members
        TypedDistribution<double>*                          branch_rate_prior;
        const TypedDagNode<Tree>*                           time_tree;
        const TypedDagNode<double>*                         root_branch_fraction;
        size_t                                              num_taxa;

    };

}

#endif
