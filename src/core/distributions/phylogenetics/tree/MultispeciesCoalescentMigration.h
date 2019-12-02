#ifndef MultispeciesCoalescentMigration_H
#define MultispeciesCoalescentMigration_H

#include "RateGenerator.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class MultispeciesCoalescentMigration : public TypedDistribution<Tree> {
        
        public:
        MultispeciesCoalescentMigration(const TypedDagNode<Tree> *st, const std::vector<Taxon> &t);
        virtual                                            ~MultispeciesCoalescentMigration(void);                                                                       //!< Virtual destructor
        
        // public member functions
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
        // public member functions
        MultispeciesCoalescentMigration*                    clone(void) const;                                                                                  //!< Create an independent clone
        void                                                setNes(TypedDagNode<RbVector<double> >* inputNes);
        void                                                setNe(TypedDagNode<double>* inputNe);
        
        protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
//        double                                              computeLnCoalescentProbability(size_t k, const std::vector<double> &t, double a, double b, size_t index, bool f);
        double                                              drawNe(size_t index);
        double                                              getNe(size_t index) const;
        double                                              getMigrationRate(size_t from, size_t to) const;

        // helper functions
        void                                                attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times);
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
//        double                                              recursivelyComputeLnProbability(const TopologyNode &n);
        void                                                resetTipAllocations(void);
        void                                                simulateTree(void);
        
        private:
        
        // parameters
        const TypedDagNode<RbVector<double> >*              Nes;
        const TypedDagNode<double >*                        Ne;
        const TypedDagNode<double >*                        delta;
        const TypedDagNode<RateGenerator>*                  Q;

        // members
        std::vector<Taxon>                                  taxa;
        const TypedDagNode<Tree>*                           species_tree;
        size_t                                              num_taxa;
        double                                              log_tree_topology_prob;
        
        std::vector< std::set< const TopologyNode* > >      individuals_per_branch;

    };
    
}

#endif
