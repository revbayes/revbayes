#ifndef EmpiricalTopologyDistribution_H
#define EmpiricalTopologyDistribution_H

#include <stddef.h>
#include <vector>

#include "Taxon.h"
#include "TraceTree.h"
#include "Tree.h"
#include "TypedDistribution.h"
#include "Clade.h"

namespace RevBayesCore {
class DagNode;
class TopologyNode;
    
    class EmpiricalTopologyDistribution : public TypedDistribution<Tree> {
        
    public:
        EmpiricalTopologyDistribution(const TraceTree& tt);
        virtual                                            ~EmpiricalTopologyDistribution(void);                                                                    //!< Virtual destructor
        
        // public member functions
        EmpiricalTopologyDistribution*                      clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                        //!< Set the current value, e.g. attach an observation (clamp)

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
    private:
        
        // helper functions
        void                                                simulateTree(void);
//        bool                                                matchesConstraints(void);
//        void                                                simulateClade(std::vector<TopologyNode*> &n);                           //!< Simulate n speciation events.

        // members
//        size_t                                              num_taxa;
//        std::vector<Taxon>                                  taxa;
//        Clade                                               outgroup;
//        bool                                                outgroup_provided;
//        bool                                                rooted;
        
        std::map<std::string, double>                       topology_probabilities;
    };
    
}

#endif
