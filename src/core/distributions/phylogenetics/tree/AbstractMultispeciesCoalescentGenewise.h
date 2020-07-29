#ifndef AbstractMultispeciesCoalescentGenewise_H
#define AbstractMultispeciesCoalescentGenewise_H

#include "RbVector.h"
#include "ModelVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    class Clade;

    class AbstractMultispeciesCoalescentGenewise : public TypedDistribution<Tree> {

    public:
        AbstractMultispeciesCoalescentGenewise(const TypedDagNode<Tree> *st, const std::vector<Taxon> &t, size_t ngt);
        virtual                                            ~AbstractMultispeciesCoalescentGenewise(void);                                                                       //!< Virtual destructor

        // public member functions
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)

        // pure virtual member functions
        virtual AbstractMultispeciesCoalescentGenewise*             clone(void) const = 0;                                                                                  //!< Create an independent clone

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
        void                                                simulateTrees(void);

        // members
        std::vector<Taxon>                                  taxa;
        const TypedDagNode<Tree>*                           species_tree;
        size_t                                              num_taxa;
        double                                              log_tree_topology_prob;
        size_t                                              num_gene_trees;
        std::vector<Tree* >*                                gene_trees;

        std::vector< std::set< const TopologyNode* > >      individuals_per_branch;

    };

}

#endif
