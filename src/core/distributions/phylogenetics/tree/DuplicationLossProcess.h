#ifndef DuplicationLossProcess_H
#define DuplicationLossProcess_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class DuplicationLossProcess : public TypedDistribution<Tree> {
        
    public:
        DuplicationLossProcess(const TypedDagNode<Tree> *it, const std::vector<Taxon> &t);
        virtual                                            ~DuplicationLossProcess(void);                                                                       //!< Virtual destructor
        
        // public member functions
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        void                                                setDuplicationRate(TypedDagNode<RbVector<double> >* d);
        void                                                setDuplicationRate(TypedDagNode<double>* d);
        void                                                setLossRate(TypedDagNode<RbVector<double> >* l);
        void                                                setLossRate(TypedDagNode<double>* l);
        void                                                setGeneSamplingProbability(TypedDagNode<RbVector<double> >* s);
        virtual void                                        setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
        // pure virtual member functions
        virtual DuplicationLossProcess*                     clone(void) const;                                                                                  //!< Create an independent clone
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
        // helper functions
        void                                                attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times);
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        double                                              computeD(double dt, double e);
        double                                              computeE(double dt, double e);
        double                                              computeLnDuplicationLossProbability(size_t k, const std::vector<double> &t, double a, double b, const TopologyNode &n, bool f);
        double                                              recursivelyComputeLnProbability(const TopologyNode &n);
        void                                                resetTipAllocations(void);
        void                                                simulateTree(void);
        
        // members
        std::vector<Taxon>                                  taxa;
        const TypedDagNode<Tree>*                           individual_tree;
        size_t                                              num_taxa;
        double                                              log_tree_topology_prob;

        const TypedDagNode<RbVector<double> >*              duplication_rates;
        const TypedDagNode<double>*                         duplication_rate;
        const TypedDagNode<RbVector<double> >*              loss_rates;
        const TypedDagNode<double>*                         loss_rate;
        const TypedDagNode<RbVector<double> >*              gene_sampling_probability;

        std::vector< std::set< const TopologyNode* > >      genes_per_branch_recent;
        std::vector< std::set< const TopologyNode* > >      genes_per_branch_ancient;
        std::vector< double >                               extinction_probs;                               // extinction probabilities at start (older/ancient end) of individual branch;

    };
    
}

#endif
