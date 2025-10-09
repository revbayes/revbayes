#ifndef UniformTopologyBranchLengthDistribution_H
#define UniformTopologyBranchLengthDistribution_H

#include <cstddef>
#include <vector>

#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"
#include "Clade.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;
    
    class UniformTopologyBranchLengthDistribution : public TypedDistribution<Tree>, public TreeChangeEventListener {
        
    public:
        UniformTopologyBranchLengthDistribution(const std::vector<Taxon> &ta, const Clade &og, TypedDistribution<double>* blp, bool rooted = false);
        UniformTopologyBranchLengthDistribution(const UniformTopologyBranchLengthDistribution &d);
        virtual                                            ~UniformTopologyBranchLengthDistribution(void);                                                                    //!< Virtual destructor

        UniformTopologyBranchLengthDistribution&            operator=(const UniformTopologyBranchLengthDistribution &d);

        // public member functions
        UniformTopologyBranchLengthDistribution*            clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        virtual void                                        fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                 //!< This node was changed in the tree
        const std::vector<Taxon>&                           getTaxa(void) const;
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                        //!< Set the current value, e.g. attach an observation (clamp)
        virtual void                                        simulateClade(std::vector<TopologyNode *> &n);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
    private:
        
        // helper functions
        void                                                simulateTree(void);
        
        // members
        TypedDistribution<double>*                          branch_length_prior;
        size_t                                              num_taxa;
        std::vector<Taxon>                                  taxa;
        double                                              log_tree_topology_prob;                                                 //!< the log-topology probability.
        Clade                                               outgroup;
        bool                                                outgroup_provided;
        bool                                                rooted;
        bool                                                dirty_topology;
        
    };
    
}

#endif

