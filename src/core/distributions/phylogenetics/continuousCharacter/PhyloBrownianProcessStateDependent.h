#ifndef PhyloBrownianProcessStateDependent_H
#define PhyloBrownianProcessStateDependent_H

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
    class PhyloBrownianProcessStateDependent : public TypedDistribution< ContinuousCharacterData > {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloBrownianProcessStateDependent(const TypedDagNode<CharacterHistoryDiscrete> *bh, size_t n_sites );
        virtual                                                            ~PhyloBrownianProcessStateDependent(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloBrownianProcessStateDependent*                         clone(void) const;                                                                      //!< Create an independent clone
        void                                                                setSigma(const TypedDagNode< double >* s);
        void                                                                setSigma(const TypedDagNode< RbVector< double > >* s);
        void                                                                setValue(ContinuousCharacterData *v, bool f=false);                                     //!< Set the current value, e.g. attach an observation (clamp)

        // non-virtual
//        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                             //!< The tree has changed and we want to know which part.
        virtual void                                                        redrawValue(void);
        double                                                              computeLnProbability(void);
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        double                                                              computeBranchTime(size_t nide_idx, double brlen);
        double                                                              computeSiteRate(size_t siteIdx);
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
        std::vector<std::vector<std::vector<double> > >                     contrasts;
        std::vector<std::vector<double> >                                   contrast_uncertainty;
        std::vector<std::vector<std::vector<double> > >                     contrast_uncertainty_per_site;
        std::vector<size_t>                                                 active_likelihood;
        
        std::vector< std::vector<bool> >                                    missing_data;

        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;
        
        bool                                                                use_missing_data;

    private:
        double                                                              computeRootState(size_t i) const;
        double                                                              computeStateDependentSigma(size_t idx) const;
        
        const TypedDagNode<CharacterHistoryDiscrete>*                       character_histories;

        const TypedDagNode< double >*                                       homogeneous_sigma;
        const TypedDagNode< RbVector< double > >*                           state_dependent_sigma;

    };
    
}


#endif

