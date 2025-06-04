#ifndef PhyloOrnsteinUhlenbeckPruning_H
#define PhyloOrnsteinUhlenbeckPruning_H

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
    class PhyloOrnsteinUhlenbeckPruning : public AbstractPhyloContinuousCharacterProcess, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloOrnsteinUhlenbeckPruning(const TypedDagNode<Tree> *t, size_t nSites );
        virtual                                                            ~PhyloOrnsteinUhlenbeckPruning(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloOrnsteinUhlenbeckPruning*                                 clone(void) const;                                                                      //!< Create an independent clone
        void                                                                setAlpha(const TypedDagNode< double >* a);
        void                                                                setAlpha(const TypedDagNode< RbVector< double > >* a);
        void                                                                setRootState(const TypedDagNode< double >* s);
        void                                                                setSigma(const TypedDagNode< double >* s);
        void                                                                setSigma(const TypedDagNode< RbVector< double > >* s);
        void                                                                setTheta(const TypedDagNode< double >* t);
        void                                                                setTheta(const TypedDagNode< RbVector< double > >* t);

        // non-virtual
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                             //!< The tree has changed and we want to know which part.
        double                                                              computeLnProbability(void);
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(const DagNode *restorer);
        void                                                                simulateRecursively(const TopologyNode& node, std::vector< ContinuousTaxonData > &t);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        double                                                              sumRootLikelihood(void);
        virtual void                                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter
        
        // the likelihoods
        std::vector<std::vector<std::vector<double> > >                     partial_likelihoods;
        std::vector<std::vector<std::vector<double> > >                     means;
        std::vector<std::vector<double> >                                   variances;
        std::vector<std::vector<std::vector<double> > >                     variances_per_site;
        std::vector<size_t>                                                 active_likelihood;

        bool                                                                use_missing_data;

        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;
        
    private:
        double                                                              computeRootState(void) const;
        double                                                              computeBranchAlpha(size_t idx) const;
        double                                                              computeBranchSigma(size_t idx) const;
        double                                                              computeBranchTheta(size_t idx) const;
        
        const TypedDagNode< double >*                                       root_state;
        const TypedDagNode< double >*                                       homogeneous_alpha;
        const TypedDagNode< double >*                                       homogeneous_sigma;
        const TypedDagNode< double >*                                       homogeneous_theta;
        const TypedDagNode< RbVector< double > >*                           heterogeneous_alpha;
        const TypedDagNode< RbVector< double > >*                           heterogeneous_sigma;
        const TypedDagNode< RbVector< double > >*                           heterogeneous_theta;

    };
    
}


#endif

