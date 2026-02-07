#ifndef PhyloMultiSampleOrnsteinUhlenbeckStateDependent_H
#define PhyloMultiSampleOrnsteinUhlenbeckStateDependent_H

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
    class PhyloMultiSampleOrnsteinUhlenbeckStateDependent : public TypedDistribution< ContinuousCharacterData > {

    public:
        enum                                                                ROOT_TREATMENT { OPTIMUM, EQUILIBRIUM, PARAMETER };
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloMultiSampleOrnsteinUhlenbeckStateDependent(const TypedDagNode<CharacterHistoryDiscrete> *bh, size_t n_sites, ROOT_TREATMENT rt, const TypedDagNode< RbVector< double > > *v, const std::vector<Taxon> &ta);
        virtual                                                            ~PhyloMultiSampleOrnsteinUhlenbeckStateDependent(void);                                             //!< Virtual destructor

        // public member functions
        // virtual
        virtual PhyloMultiSampleOrnsteinUhlenbeckStateDependent*                       clone(void) const;                                                                      //!< Create an independent clone

        void                                                                setAlpha(const TypedDagNode< double >* a);
        void                                                                setAlpha(const TypedDagNode< RbVector< double > >* a);
        void                                                                setRootValue(const TypedDagNode< double >* rvl);
        void                                                                setSigma(const TypedDagNode< double >* s);
        void                                                                setSigma(const TypedDagNode< RbVector< double > >* s);
        void                                                                setTheta(const TypedDagNode< double >* t);
        void                                                                setTheta(const TypedDagNode< RbVector< double > >* t);
        void                                                                setWithinSpeciesVariances(const TypedDagNode<RbVector<double> > *var);
        void                                                                setValue(ContinuousCharacterData *v, bool f=false);                                     //!< Set the current value, e.g. attach an observation (clamp)
        void                                                                setRootTreatment(ROOT_TREATMENT rt);
        ROOT_TREATMENT                                                      getRootTreatment() const { return root_treatment; }
        // non-virtual
        virtual void                                                        redrawValue(void);
        double                                                              computeLnProbability(void);

    protected:

        // virtual methods that may be overwritten, but then the derived class should call this methods
        double                                                              getNumberOfSamplesForSpecies(const std::string &n);
//        double                                                              computeSpeciesMean(const std::string &n, size_t i);
        double                                                              getWithinSpeciesVariance(const std::string &n);
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode &node);
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
        double                                                              computeRootValue( void ) const;
        double                                                              computeStateDependentAlpha(size_t idx) const;
        double                                                              computeStateDependentSigma(size_t idx) const;
        double                                                              computeStateDependentTheta(size_t idx) const;
        double                                                              simulateEpisode(size_t state_index, double delta_t, double ancestral_value);
        double                                                              computeEpisodeMean(double mu, size_t state_index, double time);
        double                                                              computeEpisodeScalingFactor(double log_nf, size_t state_index, double time);
        double                                                              computeEpisodeVariance(double var, size_t state_index, double time);

        ROOT_TREATMENT                                                      root_treatment;
        const TypedDagNode<CharacterHistoryDiscrete>*                       character_histories;

        const TypedDagNode< double >*                                       root_value;
        const TypedDagNode< double >*                                       homogeneous_alpha;
        const TypedDagNode< double >*                                       homogeneous_sigma;
        const TypedDagNode< double >*                                       homogeneous_theta;
        const TypedDagNode< RbVector< double > >*                           state_dependent_alpha;
        const TypedDagNode< RbVector< double > >*                           state_dependent_sigma;
        const TypedDagNode< RbVector< double > >*                           state_dependent_theta;
        const TypedDagNode< RbVector< double > >*                           within_species_variances;
        const TypedDagNode< RbVector< double > >*                           variances_of_within_species_variances;

        size_t                                                              num_species;
        size_t                                                              num_individuals;
        std::vector<size_t>                                                 num_individuals_per_species;
        std::vector<Taxon>                                                  taxa;

    };

}


#endif
