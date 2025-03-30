#ifndef PhyloBrownianProcessStateDependentTrend_H
#define PhyloBrownianProcessStateDependentTrend_H

#include "AbstractPhyloBrownianProcess.h"
#include "CharacterHistoryDiscrete.h"
#include "TreeChangeEventListener.h"

namespace RevBayesCore {
    
    /**
     * @brief A state-dependent Ornstein-Uhlenbeck process.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-23, version 1.0
     */
    class PhyloBrownianProcessStateDependentTrend : public TypedDistribution< ContinuousCharacterData > {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloBrownianProcessStateDependentTrend(const TypedDagNode<CharacterHistoryDiscrete> *bh, size_t n_sites);
        virtual                                                            ~PhyloBrownianProcessStateDependentTrend(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloBrownianProcessStateDependentTrend*                       clone(void) const;                                                                      //!< Create an independent clone
        void                                                                setRootState(const TypedDagNode< double >* s);
        void                                                                setSigma(const TypedDagNode< double >* s);
        void                                                                setSigma(const TypedDagNode< RbVector< double > >* s);
        void                                                                setTau(const TypedDagNode< double >* tau);
        void                                                                setTau(const TypedDagNode< RbVector< double > >* tau);
        void                                                                setValue(ContinuousCharacterData *v, bool f=false);                                     //!< Set the current value, e.g. attach an observation (clamp)
        // non-virtual
//        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                             //!< The tree has changed and we want to know which part.
        virtual void                                                        redrawValue(void);
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
        void                                                                simulateTipSamples( const std::vector< ContinuousTaxonData > &taxon_data );
        double                                                              sumRootLikelihood(void);
        virtual void                                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter
        
        double                                                              ln_prob;
        size_t                                                              num_nodes;
        size_t                                                              num_sites;

        // the likelihoods
        std::vector<std::vector<std::vector<double> > >                     partial_likelihoods;
        std::vector<std::vector<std::vector<double> > >                     means;
        std::vector<std::vector<double> >                                   variances;
        std::vector<size_t>                                                 active_likelihood;
        
        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;
        
    private:
        double                                                              computeRootState(void) const;
        double                                                              computeStateDependentSigma(size_t idx) const;
        double                                                              computeStateDependentTau(size_t idx) const;
        double                                                              simulateEpisode(size_t state_index, double delta_t, double ancestral_value);
        void                                                                computeEpisode(double &mu, double &variance, size_t state_index, double time);

        const TypedDagNode<CharacterHistoryDiscrete>*                       character_histories;

        const TypedDagNode< double >*                                       root_state;
        const TypedDagNode< double >*                                       homogeneous_sigma;
        const TypedDagNode< double >*                                       homogeneous_tau;
        const TypedDagNode< RbVector< double > >*                           state_dependent_sigma;
        const TypedDagNode< RbVector< double > >*                           state_dependent_tau;



    };
    
}


#endif

