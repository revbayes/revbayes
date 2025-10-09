#ifndef UniformTopologyDistribution_H
#define UniformTopologyDistribution_H

#include <cstddef>
#include <vector>

#include "Taxon.h"
#include "Tree.h"
#include "TypedDistribution.h"
#include "Clade.h"

namespace RevBayesCore {
class DagNode;
class TopologyNode;
    
    class UniformTopologyDistribution : public TypedDistribution<Tree> {
        
    public:
        UniformTopologyDistribution(const std::vector<Taxon> &ta, const Clade &og, const std::vector<Clade> &c, bool rooted = false);
		virtual                                            ~UniformTopologyDistribution(void);                                                                    //!< Virtual destructor
        
        // public member functions
        UniformTopologyDistribution*                        clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                    //!< Set the current value, e.g. attach an observation (clamp)

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // helper functions
        void                                                simulateTree(void);
        bool                                                matchesConstraints(void);
        void                                                simulateClade(std::vector<TopologyNode*> &n);                                           //!< Simulate n speciation events.

        // members
        size_t                                              num_taxa;
        std::vector<Taxon>                                  taxa;
        std::vector<Clade>                                  constraints;
        double                                              logTreeTopologyProb;                                                 //!< Topological constrains.
        Clade                                               outgroup;
        bool                                                outgroup_provided;
        bool												rooted;
    };
    
}

#endif
