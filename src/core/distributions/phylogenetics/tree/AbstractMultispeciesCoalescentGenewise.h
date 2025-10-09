#ifndef AbstractMultispeciesCoalescentGenewise_H
#define AbstractMultispeciesCoalescentGenewise_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    class Clade;

    class AbstractMultispeciesCoalescentGenewise : public TypedDistribution< RbVector<Tree> > {

    public:
        //AbstractMultispeciesCoalescentGenewise(const TypedDagNode<Tree> *st, const std::vector< std::vector<Taxon> > &t, size_t ngt);
        AbstractMultispeciesCoalescentGenewise(const TypedDagNode<Tree> *st, RbVector< RbVector<Taxon> > t, size_t ngt);
        virtual                                            ~AbstractMultispeciesCoalescentGenewise(void);                                                                       //!< Virtual destructor

        // public member functions
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        virtual void                                        setValue(RbVector<Tree>* v, bool f=false);   

        std::vector<Tree*>                                  getTrees(void) const;
        size_t                                              getNumberOfGeneTrees(void) const;

        // pure virtual member functions
        virtual AbstractMultispeciesCoalescentGenewise*             clone(void) const = 0;                                                                                  //!< Create an independent clone

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        virtual double                                      computeLnCoalescentProbability(std::vector<size_t> k, const std::vector< std::vector<double> > &t, double a, double b, size_t index, bool f) = 0;
        virtual double                                      drawNe(size_t index);

        // helper functions
        void                                                attachTimes(std::vector<Tree*> psi, std::vector< std::vector<TopologyNode *> > &tips, size_t index, const std::vector< std::vector<double> > &times);
        void                                                buildRandomBinaryTree(std::vector< std::vector<TopologyNode*> > &tips);
        double                                              recursivelyComputeLnProbability(const TopologyNode &n);
        void                                                resetTipAllocations(void);
        void                                                simulateTrees(void);

        // members
        RbVector< RbVector<Taxon> >                                         taxa;                       //!< A vector holding the vectors of taxa for each gene tree
        const TypedDagNode<Tree>*                                           species_tree;               //!< The species tree
        size_t                                                              num_species;                //!< The number of tips in the species tree
        std::vector<size_t>                                                 num_taxa;                   //!< A vector holding the number of tips for each gene tree
        double                                                              log_tree_topology_prob;     //!< Combinatorial topology prob for species tree
        size_t                                                              num_gene_trees;             //!< Number of genes/gene trees
        std::vector<Tree*>                                                  gene_trees;                 //!< Vector of gene trees

        std::vector< std::vector< std::set< const TopologyNode* > > >       individuals_per_branch_genewise;     //!< A vector holding the vectors that contain the individuals per branch for each gene tree

    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const AbstractMultispeciesCoalescentGenewise& x);

}

#endif
