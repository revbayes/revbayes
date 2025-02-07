#ifndef BranchRateTreeDistribution_H
#define BranchRateTreeDistribution_H

#include <cstddef>
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
        RbBitSet                                            collectTreeSample(const Tree& t, const TopologyNode& n, RbBitSet& in, std::map<RbBitSet, double>& bl);
        RbBitSet                                            collectSplits(const TopologyNode& n, RbBitSet& in, std::vector<RbBitSet>& s) const;
        void                                                simulateTree(void);
        void                                                simulateClade(std::vector<TopologyNode*> &n);                           //!< Simulate n speciation events.
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // members
        TypedDistribution<double>*                          branch_rate_prior;
        const TypedDagNode<Tree>*                           time_tree;
        const TypedDagNode<double>*                         root_branch_fraction;
        size_t                                              num_taxa;
        
        // variables catching some of the probability computation
        bool                                                touched_time_tree;
        bool                                                touched_branch_length_tree;
        bool                                                was_touched_time_tree;
        bool                                                was_touched_branch_length_tree;
        std::string                                         newick_time_tree;
        std::string                                         newick_branch_length_tree;
        Tree*                                               time_tree_unrooted;
        std::vector<RbBitSet>                               splits;
        std::map<RbBitSet, double>                          split_to_branch_lengths;
        std::string                                         stored_newick_time_tree;
        std::string                                         stored_newick_branch_length_tree;
        Tree*                                               stored_time_tree_unrooted;
        std::vector<RbBitSet>                               stored_splits;
        std::map<RbBitSet, double>                          stored_split_to_branch_lengths;


    };

}

#endif
