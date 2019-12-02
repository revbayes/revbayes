#ifndef PhylogeneticIndependentContrastsFunction_H
#define PhylogeneticIndependentContrastsFunction_H

#include "ContinuousCharacterData.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class PhylogeneticIndependentContrastsFunction : public TypedFunction< RbVector<double> > {
        
    public:
        PhylogeneticIndependentContrastsFunction(const TypedDagNode<Tree> *t, const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<long>* s, bool n);
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
        const TypedDagNode<long>*                               site;
        bool                                                    normalized;
        
        std::vector<double>                                     contrasts;
        std::vector<double>                                     standard_deviations;
        std::vector<double>                                     node_means;
        std::vector<double>                                     contrast_uncertainty;
    };
    
}

#endif

