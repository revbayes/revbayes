#ifndef PhylogeneticIndependentContrastsMultiSampleFunction_H
#define PhylogeneticIndependentContrastsMultiSampleFunction_H

#include "ContinuousCharacterData.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class PhylogeneticIndependentContrastsMultiSampleFunction : public TypedFunction< RbVector<double> > {
        
    public:
        PhylogeneticIndependentContrastsMultiSampleFunction(const TypedDagNode<Tree> *t, const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<long>* s, const std::vector<Taxon> &ta, bool n);
        virtual                                                ~PhylogeneticIndependentContrastsMultiSampleFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        PhylogeneticIndependentContrastsMultiSampleFunction*                clone(void) const;                                                                  //!< Create an independent clone
        void                                                                update(void);
        
    protected:
        double                                                              computeMeanForSpecies(const std::string &n, size_t i);
        double                                                              getNumberOfSamplesForSpecies(const std::string &n);
        double                                                              getWithinSpeciesVariance(const std::string &n, size_t i);
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                resetContrasts(void);
        
    private:
        
        // members
        const TypedDagNode<Tree>*                                           tau;
        const TypedDagNode<ContinuousCharacterData>*                        data;
        const TypedDagNode<long>*                                           site;
        bool                                                                normalized;
        std::vector<Taxon>                                                  taxa;

        std::vector<double>                                                 contrasts;
        std::vector<double>                                                 standard_deviations;
        std::vector<double>                                                 node_means;
        std::vector<double>                                                 contrast_uncertainty;
    };
    
}

#endif


