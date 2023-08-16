#ifndef PhyloMultivariateBrownianProcessMultiSampleREML_H
#define PhyloMultivariateBrownianProcessMultiSampleREML_H

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
     * @since 2019-04-29, version 1.0.11
     */
    class PhyloMultivariateBrownianProcessMultiSampleREML : public AbstractPhyloBrownianProcess, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
//        PhyloMultivariateBrownianProcessMultiSampleREML(const TypedDagNode<Tree> *tr, const TypedDagNode<MatrixReal> *c, const TypedDagNode< RbVector< double > > *v, const TypedDagNode< RbVector< double > > *v2, const std::vector<Taxon> &ta, size_t ns );
        PhyloMultivariateBrownianProcessMultiSampleREML(const TypedDagNode<Tree> *tr, const TypedDagNode<MatrixReal> *c, const TypedDagNode< RbVector< RbVector< double > > > *v, const std::vector<Taxon> &ta, size_t ns );
        virtual                                                            ~PhyloMultivariateBrownianProcessMultiSampleREML(void);                                                                      //!< Virtual destructor
        
        // public member functions
        virtual PhyloMultivariateBrownianProcessMultiSampleREML*            clone(void) const;                                                                                              //!< Create an independent clone
        double                                                              computeLnProbability(void);
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                //!< The tree has changed and we want to know which part.
        void                                                                redrawValue(void);
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        double                                                              computeMeanForSpecies(const std::string &n, size_t i);
        double                                                              getNumberOfSamplesForSpecies(const std::string &n);
        double                                                              getWithinSpeciesVariance(const std::string &n, size_t site);
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(const DagNode *restorer);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        void                                                                simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa);
        virtual void                                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                                //!< Swap a parameter
        
        const TypedDagNode< RbVector< RbVector< double > > >*               within_species_variances;

        
        // the likelihoods
        std::vector<std::vector<double> >                                   partial_branch_likelihoods;
        std::vector<std::vector<std::vector<double> > >                     contrasts;
        std::vector<std::vector<std::vector<double> > >                     contrast_uncertainty;
        std::vector<size_t>                                                 active_likelihood;
        
        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;
        std::vector<size_t>                                                 site_indices;
        
        
        std::vector<Taxon>                                                  taxa;
        const TypedDagNode< MatrixReal >*                                   rate_matrix;
        size_t                                                              active_matrix;
        std::vector<MatrixReal>                                             precision_matrices;
        
        
    };
    
}


#endif
