#ifndef PhyloBrownianProcessMultiSampleREML_H
#define PhyloBrownianProcessMultiSampleREML_H

#include "AbstractPhyloBrownianProcess.h"
#include "IndexedCache.h"
#include "TreeChangeEventListener.h"

namespace RevBayesCore {
    
    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-23, version 1.0
     */
    class PhyloBrownianProcessMultiSampleREML : public AbstractPhyloBrownianProcess, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloBrownianProcessMultiSampleREML(const TypedDagNode<Tree> *tr, const TypedDagNode< RbVector< double > > *v, const std::vector<Taxon> &ta, size_t ns );
        virtual                                                            ~PhyloBrownianProcessMultiSampleREML(void);                                                                      //!< Virtual destructor
        
        // public member functions
        virtual PhyloBrownianProcessMultiSampleREML*                        clone(void) const;                                                                                              //!< Create an independent clone
        double                                                              computeLnProbability(void);
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                //!< The tree has changed and we want to know which part.
        void                                                                redrawValue(void);
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        double                                                              computeMeanForSpecies(const std::string &n, size_t i);
        double                                                              getNumberOfSamplesForSpecies(const std::string &n);
        double                                                              getWithinSpeciesVariance(const std::string &n);
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        void                                                                invalidateBranchAndAncestors(const TopologyNode& n);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(const DagNode *restorer);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        virtual void                                                        snapshotSpecialization(void);
        double                                                              sumRootLikelihood(void);
        virtual void                                                        invalidateSpecialization(const DagNode *toucher, bool touchAll);
        
        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                                //!< Swap a parameter
        
        const TypedDagNode< RbVector< double > >*                           within_species_variances;

        
        struct NodeCache
        {
            std::vector<double>                                             partial_likelihoods;
            std::vector<double>                                             means;
            std::vector<double>                                             variances;
            std::vector<bool>                                               missing_data;
        };

        IndexedCache<NodeCache>                                             node_likelihoods;
        std::vector<size_t>                                                 site_indices;

        
        std::vector<Taxon>                                                  taxa;
        std::map<string,size_t>                                             sample_to_species_index;

    };
    
}


#endif
