//
//  MixedSamplingBirthDeathProcess.h
//  revbayes-proj
//
//  Created by Michael Landis on 11/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef MixedSamplingBirthDeathProcess_h
#define MixedSamplingBirthDeathProcess_h


#include "Clade.h"
#include "RbVector.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * @file
     * This file contains the declaration of the random variable class for constant rate birth-death process.
     *
     * @brief Declaration of the constant rate Birth-Death process class.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-17, version 1.0
     *
     */
    class MixedSamplingBirthDeathProcess : public TypedDistribution<Tree>, TreeChangeEventListener {
        
    public:
        MixedSamplingBirthDeathProcess(TypedDistribution<Tree>* base_dist_1, TypedDistribution<Tree>* base_dist_2, const TypedDagNode<Tree>* bb);
        MixedSamplingBirthDeathProcess(const MixedSamplingBirthDeathProcess &d);
        
        virtual ~MixedSamplingBirthDeathProcess(void);
        // pure virtual member functions
        virtual MixedSamplingBirthDeathProcess*        clone(void) const;                                                                                  //!< Create an independent clone
        
        
        // public member functions you may want to override
        double                                              computeLnProbability(void);                                                                         //!< Compute the log-transformed probability of the current value.
        void                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        virtual void                                        redrawValue(void);                                                                                  //!< Draw a new random value from the distribution
//        void                                                setBackbone( const TypedDagNode<Tree> *backbone_one=NULL, const TypedDagNode<RbVector<Tree> > *backbone_many=NULL);
        virtual void                                        setStochasticNode(StochasticNode<Tree> *n);                                                         //!< Set the stochastic node holding this distribution
        virtual void                                        setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
        
    protected:
        
        void                                                initializeBitSets();
        virtual void                                        keepSpecialization(DagNode* affecter);
        virtual void                                        restoreSpecialization(DagNode *restorer);
        virtual void                                        touchSpecialization(DagNode *toucher, bool touchAll);
        
        // Parameter management functions. You need to override both if you have additional parameters
        virtual void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        
        
        // helper functions
        bool                                                matchesBackbone(void);
//        bool                                                matchesConstraints(void);
        RbBitSet                                            recursivelyAddBackboneConstraints(const TopologyNode& node, size_t backbone_idx);
        void                                                recursivelyFlagNodesDirty(const TopologyNode& n);
        RbBitSet                                            recursivelyUpdateClades(const TopologyNode& node);
        Tree*                                               simulateTree(void);
        
        
        // members
        std::vector<std::vector<RbBitSet> >                 active_backbone_clades;
        std::vector<RbBitSet>                               active_clades;
        std::vector<std::vector<RbBitSet> >                 backbone_constraints;
        std::vector<RbBitSet>                               backbone_mask;
        
        const TypedDagNode<Tree>*                           backbone_topology;
//        const TypedDagNode<RbVector<Tree> >*                backbone_topologies;
        
        TypedDistribution<Tree>*                            base_distribution_1;
        TypedDistribution<Tree>*                            base_distribution_2;
        std::vector<bool>                                   dirty_nodes;
//        std::vector<Clade>                                  monophyly_constraints;
        std::vector<std::vector<RbBitSet> >                 stored_backbone_clades;
        std::vector<RbBitSet>                               stored_clades;
        size_t                                              num_backbones;
//        bool                                                use_multiple_backbones;
        
    };
    
}

#endif /* MixedSamplingBirthDeathProcess_h */
