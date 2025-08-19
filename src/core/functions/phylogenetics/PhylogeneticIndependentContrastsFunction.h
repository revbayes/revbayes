#ifndef PhylogeneticIndependentContrastsFunction_H
#define PhylogeneticIndependentContrastsFunction_H

#include <cstddef>
#include <vector>

#include "RbVector.h"
#include "TypedFunction.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class ContinuousCharacterData;
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class PhylogeneticIndependentContrastsFunction : public TypedFunction< RbVector<double> > {
        
    public:
        PhylogeneticIndependentContrastsFunction(const TypedDagNode<Tree> *t, const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t>* s, bool n);
        virtual                                                ~PhylogeneticIndependentContrastsFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        PhylogeneticIndependentContrastsFunction*               clone(void) const;                                                                  //!< Create an independent clone
        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        void                                                    recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                    resetContrasts(void);
        
    private:
        
        // members
        const TypedDagNode<Tree>*                               tau;
        const TypedDagNode<ContinuousCharacterData>*            data;
        const TypedDagNode<std::int64_t>*                               site;
        bool                                                    normalized;
        
        std::vector<double>                                     contrasts;
        std::vector<double>                                     standard_deviations;
        std::vector<double>                                     node_means;
        std::vector<double>                                     contrast_uncertainty;
    };
    
}

#endif

