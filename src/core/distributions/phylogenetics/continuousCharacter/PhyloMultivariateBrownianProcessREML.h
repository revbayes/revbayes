#ifndef PhyloMultivariateBrownianProcessREML_H
#define PhyloMultivariateBrownianProcessREML_H

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
    class PhyloMultivariateBrownianProcessREML : public AbstractPhyloBrownianProcess, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloMultivariateBrownianProcessREML(const TypedDagNode<Tree> *t, const TypedDagNode<MatrixReal> *c, size_t nSites );
        virtual                                                            ~PhyloMultivariateBrownianProcessREML(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloMultivariateBrownianProcessREML*                       clone(void) const;                                                                      //!< Create an independent clone
        
        // non-virtual
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                             //!< The tree has changed and we want to know which part.
        double                                                              computeLnProbability(void);
        std::vector<std::vector<double> >                                   getContrasts(void);
        std::vector<double>                                                 getContrastStDevs(void){ return independent_contrasts_sds;}
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        void                                                                invalidateBranchAndAncestors(const TopologyNode& n);
        void                                                                invalidateInternalNodes(void);
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursiveComputeContrasts( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(const DagNode *restorer);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        virtual void                                                        simulateRecursively(const TopologyNode& node, std::vector< ContinuousTaxonData > &t);
        virtual void                                                        snapshotSpecialization(void);
        virtual void                                                        invalidateSpecialization(const DagNode *toucher, bool touchAll);

        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter

        struct NodeCache
        {
            double                                                          partial_likelihood = 0.0;
            std::vector<double>                                             contrasts;
            double                                                          contrast_uncertainty = 0.0;
        };

        IndexedCache<NodeCache>                                             node_likelihoods;
        
        std::vector<std::vector<double> >                                   independent_contrasts;
        std::vector<double>                                                 independent_contrasts_sds;

    private:
        
        const TypedDagNode< MatrixReal >*                                   rate_matrix;
        size_t                                                              active_matrix;
        std::vector<MatrixReal>                                             precision_matrices;
        
    };
    
}


#endif
