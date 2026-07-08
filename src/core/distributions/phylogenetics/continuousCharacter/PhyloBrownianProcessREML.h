#ifndef PhyloBrownianProcessREML_H
#define PhyloBrownianProcessREML_H

#include "AbstractPhyloBrownianProcess.h"
#include "SnapshotCache.h"
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
    class PhyloBrownianProcessREML : public AbstractPhyloBrownianProcess, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloBrownianProcessREML(const TypedDagNode<Tree> *t, size_t nSites );
        virtual                                                            ~PhyloBrownianProcessREML(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloBrownianProcessREML*                                   clone(void) const;                                                                      //!< Create an independent clone
        
        // non-virtual
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                             //!< The tree has changed and we want to know which part.
        double                                                              computeLnProbability(void);
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                        keepSpecialization(void);
        void                                                                invalidateBranchAndAncestors(const TopologyNode& n);
        void                                                                invalidateInternalNodes(void);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(void);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        virtual void                                                        snapshotSpecialization(void);
        double                                                              sumRootLikelihood(void);
        virtual void                                                        invalidateSpecialization(const DagNode *toucher, bool fullyInvalidateSelf);

        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter

        struct NodeCache
        {
            std::vector<double>                                             partial_likelihoods;
            std::vector<double>                                             means;
            std::vector<double>                                             variances_per_site;
            std::vector<bool>                                               missing_data;
            double                                                          variance = 0.0;
        };

        IndexedSnapshotCache<NodeCache>                                             node_likelihoods;
        
    private:
                
    };
    
}


#endif
