#ifndef TopologyConstrainedTreeDistribution_H
#define TopologyConstrainedTreeDistribution_H

#include <cstdint>

#include "Clade.h"
#include "RbVector.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     */
    class TopologyConstrainedTreeDistribution : public TypedDistribution<Tree>, TreeChangeEventListener {
        
    public:
        TopologyConstrainedTreeDistribution(TypedDistribution<Tree>* base_dist, const std::vector<Clade> &c, Tree *t, std::int64_t age_check_precision);
        TopologyConstrainedTreeDistribution(const TopologyConstrainedTreeDistribution &d);
        
        virtual ~TopologyConstrainedTreeDistribution(void);
        
        TopologyConstrainedTreeDistribution&                operator=(const TopologyConstrainedTreeDistribution &d);

        
        // pure virtual member functions
        virtual TopologyConstrainedTreeDistribution*        clone(void) const;                                                                                  //!< Create an independent clone
        
        
        // public member functions you may want to override
        double                                              computeLnProbability(void);                                                                         //!< Compute the log-transformed probability of the current value.
        void                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                    //!< The tree has changed and we want to know which part.
        virtual void                                        redrawValue(SimulationCondition c);                                                                 //!< Draw a new random value from the distribution
        virtual void                                        redrawValue(void);                                                                                  //!< Draw a new random value from the distribution
        
        void                                                setBackbone( const TypedDagNode<Tree> *backbone_one=NULL, const TypedDagNode<RbVector<Tree> > *backbone_many=NULL);
        virtual void                                        setStochasticNode(StochasticNode<Tree> *n);                                                         //!< Set the stochastic node holding this distribution
        virtual void                                        setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        virtual bool                                        allowsSA(void) { return base_distribution->allowsSA(); }                                            //!< Checks if distribution is compatible with sampled ancestors
        
    protected:
        
        void                                                initializeBitSets();
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                  //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // Parameter management functions. You need to override both if you have additional parameters
        virtual void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        
        
        // helper functions
        bool                                                matchesBackbone(void);
        bool                                                matchesConstraints(void);
        RbBitSet                                            recursivelyAddBackboneConstraints(const TopologyNode& node, size_t backbone_idx);
        void                                                recursivelyFlagNodesDirty(const TopologyNode& n);
        RbBitSet                                            recursivelyUpdateClades(const TopologyNode& node);
        Tree*                                               simulateRootedTree(bool alwaysReturn);
        Tree*                                               simulateUnrootedTree(void);


        // members
        std::vector<std::vector<RbBitSet> >                 active_backbone_clades;
        std::vector<RbBitSet>                               active_clades;
        std::vector<std::vector<RbBitSet> >                 backbone_constraints;
        std::vector<RbBitSet>                               backbone_mask;
        
        const TypedDagNode<Tree>*                           backbone_topology;
        const TypedDagNode<RbVector<Tree> >*                backbone_topologies;
        
        TypedDistribution<Tree>*                            base_distribution;
        std::vector<bool>                                   dirty_nodes;
        std::vector<Clade>                                  monophyly_constraints;
        std::vector<std::vector<RbBitSet> >                 stored_backbone_clades;
        std::vector<RbBitSet>                               stored_clades;
        size_t                                              num_backbones;
        bool                                                use_multiple_backbones;
        Tree*                                               starting_tree;
        
        bool                                                rooting_known;
        bool                                                is_rooted;

    };
    
}

#endif
