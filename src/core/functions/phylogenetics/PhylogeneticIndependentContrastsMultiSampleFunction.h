#ifndef PhylogeneticIndependentContrastsMultiSampleFunction_H
#define PhylogeneticIndependentContrastsMultiSampleFunction_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "RbVector.h"
#include "TypedFunction.h"
#include "Taxon.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class ContinuousCharacterData;
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class PhylogeneticIndependentContrastsMultiSampleFunction : public TypedFunction< RbVector<double> > {
        
    public:
        PhylogeneticIndependentContrastsMultiSampleFunction(const TypedDagNode<Tree> *t, const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t>* s, const std::vector<Taxon> &ta, bool n);
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
        const TypedDagNode<std::int64_t>*                                           site;
        bool                                                                normalized;
        std::vector<Taxon>                                                  taxa;

        std::vector<double>                                                 contrasts;
        std::vector<double>                                                 standard_deviations;
        std::vector<double>                                                 node_means;
        std::vector<double>                                                 contrast_uncertainty;
    };
    
}

#endif


