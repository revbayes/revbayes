#ifndef PhyloMultivariateBrownianProcessREML_H
#define PhyloMultivariateBrownianProcessREML_H

#include "AbstractPhyloBrownianProcess.h"
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
        virtual void                                                        keepSpecialization(DagNode* affecter);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursiveComputeContrasts( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(DagNode *restorer);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        virtual void                                                        simulateRecursively(const TopologyNode& node, std::vector< ContinuousTaxonData > &t);
        virtual void                                                        touchSpecialization(DagNode *toucher, bool touchAll);

        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter

        // the likelihoods
        std::vector<std::vector<double> >                                   partial_likelihoods;
        std::vector<std::vector<std::vector<double> > >                     contrasts;
        std::vector<std::vector<double> >                                   contrast_uncertainty;
        std::vector<size_t>                                                 active_likelihood;
        
        std::vector<std::vector<double> >                                   independent_contrasts;
        std::vector<double>                                                 independent_contrasts_sds;
        
        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;

    private:
        
        const TypedDagNode< MatrixReal >*                                   rate_matrix;
        size_t                                                              active_matrix;
        std::vector<MatrixReal>                                             precision_matrices;
        
    };
    
}


#endif
