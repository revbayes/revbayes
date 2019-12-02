#ifndef AbstractMultispeciesCoalescent_H
#define AbstractMultispeciesCoalescent_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class AbstractMultispeciesCoalescent : public TypedDistribution<Tree> {
        
    public:
        AbstractMultispeciesCoalescent(const TypedDagNode<Tree> *st, const std::vector<Taxon> &t);
        virtual                                            ~AbstractMultispeciesCoalescent(void);                                                                       //!< Virtual destructor
        
        // public member functions
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)

        // pure virtual member functions
        virtual AbstractMultispeciesCoalescent*             clone(void) const = 0;                                                                                  //!< Create an independent clone

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        virtual double                                      computeLnCoalescentProbability(size_t k, const std::vector<double> &t, double a, double b, size_t index, bool f) = 0;
        virtual double                                      drawNe(size_t index);

        // helper functions
        void                                                attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times);
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        double                                              recursivelyComputeLnProbability(const TopologyNode &n);
        void                                                resetTipAllocations(void);
        void                                                simulateTree(void);
        
        // members
        std::vector<Taxon>                                  taxa;
        const TypedDagNode<Tree>*                           species_tree;
        size_t                                              num_taxa;
        double                                              log_tree_topology_prob;
        
        std::vector< std::set< const TopologyNode* > >      individuals_per_branch;

    };
    
}

#endif
