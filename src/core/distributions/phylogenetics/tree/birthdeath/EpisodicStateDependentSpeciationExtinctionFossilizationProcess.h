#ifndef EpisodicStateDependentSpeciationExtinctionFossilizationProcess_H
#define EpisodicStateDependentSpeciationExtinctionFossilizationProcess_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TreeDiscreteCharacterData.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "RateMatrix.h"
#include "RateMatrix_JC.h"
#include "Simplex.h"
#include "SSE_ODE.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDagNode.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlString.h"


#include <vector>

namespace RevBayesCore {
    
    class Clade;
    
    /**
     * @file
     * This file contains the declaration of the random variable class for the character-dependent
     * cladogenetic birth-death process: ClaSSE as described in Golberg & Igic 2012
     *
     * Will Freyman 6/22/16
     *
     */
    class EpisodicStateDependentSpeciationExtinctionFossilizationProcess : public TypedDistribution<Tree>, public TreeChangeEventListener, public MemberObject< RbVector<std::int64_t> >, public MemberObject< RbVector<double> > {
        
    public:
        EpisodicStateDependentSpeciationExtinctionFossilizationProcess(const TypedDagNode<double> *root,
                                                                       const TypedDagNode<Simplex>* p,
                                                                       const std::string &cdt,
                                                                       bool uo,
                                                                       size_t min_num_lineages,
                                                                       size_t max_num_lineages,
                                                                       size_t exact_num_lineages,
                                                                       double max_t,
                                                                       bool prune,
                                                                       bool condition_on_tip_states,
                                                                       bool condition_on_num_tips,
                                                                       bool condition_on_tree,
                                                                       std::int64_t age_check_precision);
        
        // pure virtual member functions
        virtual EpisodicStateDependentSpeciationExtinctionFossilizationProcess*              clone(void) const;
        virtual                                                         ~EpisodicStateDependentSpeciationExtinctionFossilizationProcess(void);                                                              //!< Virtual destructor

        double                                                          computeLnProbability(void);
        void                                                            fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        const AbstractHomologousDiscreteCharacterData&                  getCharacterData() const;
        double                                                          getOriginAge(void) const;
        std::vector<double>                                             getAverageExtinctionRatePerBranch(void) const;
        std::vector<double>                                             getAverageSpeciationRatePerBranch(void) const;
        std::vector<std::int64_t>                                       getNumberOfShiftEventsPerBranch(void) const;
        std::vector<double>                                             getTimeInStates(void) const;
        double                                                          getRootAge(void) const;
        virtual void                                                    redrawValue(void);
        void                                                            setCladogenesisMatrix(const TypedDagNode< CladogeneticSpeciationRateMatrix > *r);
//        void                                                            setCladogenesisMatrix(const TypedDagNode< CladogeneticSpeciationRateMatrix > *r);
        void                                                            setExtinctionRates(const TypedDagNode< RbVector<double> > *r);
        void                                                            setExtinctionRates(const TypedDagNode< RbVector< RbVector<double> > > *r, const TypedDagNode<RbVector<double> >* t);
        void                                                            setFossilizationRates(const TypedDagNode< RbVector<double> > *r);
        void                                                            setFossilizationRates(const TypedDagNode< RbVector<RbVector<double> > > *r, const TypedDagNode<RbVector<double> >* t);
        void                                                            setMassExtinctionSurvivalProbabilities(const TypedDagNode<RbVector<RbVector<double> > > *p, const TypedDagNode<RbVector<double> >* t);
        void                                                            setSampleCharacterHistory(bool sample_history);                                                     //!< Set whether or not we are sampling the character history along branches.
        void                                                            setSamplingFraction(const TypedDagNode< double > *r);
        void                                                            setSamplingFraction(const TypedDagNode< RbVector<double> > *r);
        void                                                            setSpeciationRates(const TypedDagNode< RbVector<double> > *r);
        void                                                            setSpeciationRates(const TypedDagNode< RbVector< RbVector<double> > > *r, const TypedDagNode<RbVector<double> >* t);
        void                                                            setTransitionRate(const TypedDagNode<double> *r);
        void                                                            setTransitionRate(const TypedDagNode< RbVector<double> > *r, const TypedDagNode<RbVector<double> >* t);
        void                                                            setTransitionRateMatrix(const TypedDagNode< RateGenerator > *m);
        void                                                            setTransitionRateMatrix(const TypedDagNode< RbVector< RateGenerator > > *m, const TypedDagNode<RbVector<double> >* t);
        void                                                            setNumberOfTimeSlices(double n);                                                                    //!< Set the number of time slices for the numerical ODE.
        virtual void                                                    setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
        void                                                            drawJointConditionalAncestralStates(std::vector<size_t>& startStates, std::vector<size_t>& endStates);
        void                                                            recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<size_t>& startStates, std::vector<size_t>& endStates);
        void                                                            drawStochasticCharacterMap(std::vector<std::string>& character_histories, bool set_amb_char_data = false, bool use_simmap_default=true);
        bool                                                            recursivelyDrawStochasticCharacterMap(const TopologyNode &node, size_t start_state, std::vector<std::string>& character_histories, bool set_amb_char_data, bool use_simmap_default);
        void                                                            numericallyIntegrateProcess(std::vector< double > &likelihoods, double begin_age, double end_age, bool use_backward, bool extinction_only) const; //!< Wrapper function for the ODE time stepper function.
        void                                                            resizeVectors(size_t num_nodes);

    protected:
        
        double                                                          getEventRate(double a) const;
        const RateGenerator&                                            getEventRateMatrix(double a) const;
        std::vector<double>                                             getRootFrequencies(void) const;

        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                    getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                  //!< get affected nodes
        virtual void                                                    keepSpecialization(const DagNode* affecter);
        virtual void                                                    restoreSpecialization(const DagNode *restorer);
        virtual void                                                    touchSpecialization(const DagNode *toucher, bool touchAll);
        
        double                                                          lnProbTreeShape(void) const;
        
        // Parameter management functions. You need to override both if you have additional parameters
        virtual void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        void                                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;
        void                                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<std::int64_t> &rv) const;     //!< Map the member methods to internal function calls
        RevLanguage::RevPtr<RevLanguage::RevVariable>                   executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found);
        
        // helper functions
        void                                                            addTimesToGlobalTimeline(std::set<double> &event_times, const RbVector<double>& par_times) const;        //!< Adds timeline for parameter to set that we will use for global timeline
        void                                                            buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        void                                                            checkVectorSizes(const TypedDagNode<RbVector<double> >* v1, const TypedDagNode<RbVector<double> >* v2, int v1_minus_v2, const std::string& param_name, bool is_rate) const;
        std::vector<double>                                             pExtinction(double start, double end) const;                                                        //!< Compute the probability of extinction of the process (without incomplete taxon sampling).
        virtual double                                                  pSurvival(double start, double end) const;                                                          //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double                                                          pSurvival(double start, double end, bool speciation) const;                                                          //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        void                                                            recursivelyFlagNodeDirty(const TopologyNode& n);
        bool                                                            simulateTree(size_t attempts = 0);
        bool                                                            simulateTreeConditionedOnTips(size_t attempts = 0);
        std::vector<double>                                             calculateExtinctionRatePerState(double a) const;
        std::vector<double>                                             calculateTotalAnageneticRatePerState(double a) const;
        std::vector<double>                                             calculateTotalSpeciationRatePerState(double a) const;
        double                                                          computeEpochBegin(size_t i) const;
        size_t                                                          computeEpochIndex(double a) const;
        double                                                          computeEpochEnd(size_t i) const;
        void                                                            computeNodeProbability(const TopologyNode &n, size_t nIdx) const;
        double                                                          computeRootLikelihood() const;
        const RbVector<double>&                                         computeExtinctionRateAtTime(double a) const;
        const RbVector<double>&                                         computeSurvivalProbabilitiesAtTime(double a) const;
        const RbVector<double>&                                         computeFossilizationRateAtTime(double a) const;
        const RbVector<double>&                                         computeSpeciationRateAtTime(double a) const;
        void                                                            expandNonGlobalProbabilityParameterVector(std::vector< RbVector<double> > &par, const std::vector<double> &par_times, const RbVector<double> &d) const; //!< Updates vector par such that it matches the global timeline
        void                                                            expandNonGlobalRateParameterVector(std::vector<RbVector<double> > &par, const std::vector<double> &par_times) const; //!< Updates vector par such that it matches the global timeline
        void                                                            expandNonGlobalRateParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const;   //!< Updates vector par such that it matches the global timeline
        void                                                            expandNonGlobalRateParameterVector(std::vector<size_t> &par, const std::vector<double> &par_times) const;   //!< Updates vector par such that it matches the global timeline
        size_t                                                          findIndex(double t) const;                                                                                  //!< Find the index so that times[index-1] < t < times[index]
        size_t                                                          findIndex(double t, const std::vector<double>& timeline) const;
        bool                                                            isEpisodicModel(void) const;                                                                                //!< Checks if we have a constant-rate process
        void                                                            prepareTimeline(void) const;
        void                                                            sortGlobalTimesAndVectorParameter(void) const;                                                              //!< Sorts times to run from 0->inf, and orders ALL vector parameters to match
        void                                                            sortNonGlobalTimesAndParameters(std::vector<RbVector<double> >& par, std::vector<double>& times) const;     //!< Sorts times to run from 0->inf, and orders par to match
        void                                                            sortNonGlobalTimesAndParameters(std::vector<double>& par, std::vector<double>& times) const;                //!< Sorts times to run from 0->inf, and orders par to match
        void                                                            sortNonGlobalTimesAndParameters(std::vector<size_t>& par, std::vector<double>& times) const;           //!< Sorts times to run from 0->inf, and orders par to match

        // members
        std::string                                                     condition;                                  //!< The condition of the process (none/survival/#taxa).
        double                                                          dt;                                         //!< The size of the time slices used by the ODE for numerical integration.
        std::vector<bool>                                               active_likelihood;
        mutable std::vector<bool>                                       changed_nodes;
        mutable std::vector<bool>                                       dirty_nodes;
        mutable std::vector<std::vector<std::vector<double> > >         node_partial_likelihoods;
        mutable std::map<size_t, std::vector<std::vector<double> > >    branch_partial_likelihoods;
        mutable std::vector<std::vector<double> >                       extinction_probabilities;
        size_t                                                          num_states;
        mutable std::vector<std::vector<double> >                       scaling_factors;
        bool                                                            use_cladogenetic_events;                    //!< do we use the speciation rates from the cladogenetic event map?
        mutable bool                                                    use_episodic_model;                         //!< do we use the speciation rates from the cladogenetic event map?
        bool                                                            use_origin;
        bool                                                            sample_character_history;                   //!< are we sampling the character history along branches?
        std::vector<double>                                             average_speciation;
        std::vector<double>                                             average_extinction;
        std::vector<std::int64_t>                                       num_shift_events;
        std::vector<double>                                             time_in_states;
        std::string                                                     simmap;
        
        // parameters
        const TypedDagNode< CladogeneticSpeciationRateMatrix >*         cladogenesis_matrix;
        const TypedDagNode<double>*                                     process_age;                                //!< Time since the origin.
        const TypedDagNode<RbVector<double> >*                          mu_const;
        const TypedDagNode<RbVector<RbVector<double> > >*               mu_var;
        const TypedDagNode<RbVector<double> >*                          lambda_const;
        const TypedDagNode<RbVector<RbVector<double> > >*               lambda_var;
        const TypedDagNode<RbVector<double> >*                          phi_const;
        const TypedDagNode<RbVector<RbVector<double> > >*               phi_var;
        const TypedDagNode<double>*                                     eta_const;
        const TypedDagNode<RbVector<double> >*                          eta_var;
        const TypedDagNode<RateGenerator>*                              Q_const;
        const TypedDagNode<RbVector<RateGenerator> >*                   Q_var;
        const TypedDagNode<RbVector<RbVector<double> > >*               survival_probs;
        const TypedDagNode<RbVector<double> >*                          epoch_times_lambda;
        const TypedDagNode<RbVector<double> >*                          epoch_times_mu;
        const TypedDagNode<RbVector<double> >*                          epoch_times_phi;
        const TypedDagNode<RbVector<double> >*                          epoch_times_gamma;
        const TypedDagNode<RbVector<double> >*                          epoch_times_Q;
        const TypedDagNode<RbVector<double> >*                          epoch_times_eta;
        const TypedDagNode<Simplex >*                                   pi;                                         //!< The root frequencies (probabilities of the root states).
        const TypedDagNode<double>*                                     rho;                                        //!< Sampling probability of each species.
        const TypedDagNode<RbVector<double> >*                          rho_per_state;                              //!< Sampling probability of each species.

        mutable std::vector<RbVector<double> >                          lambda;
        mutable std::vector<RbVector<double> >                          mu;
        mutable std::vector<RbVector<double> >                          phi;
        mutable std::vector<RbVector<double> >                          gamma;
        mutable std::vector<double>                                     eta;
        mutable RbVector<RateGenerator>                                 Q;
        mutable std::vector<size_t>                                     Q_indices;

        mutable std::vector<double>                                     global_timeline;                            //!< The times of the instantaneous events and rate shifts.

        
        RateMatrix_JC                                                   Q_default;
        size_t                                                          min_num_lineages;
        size_t                                                          max_num_lineages;
        size_t                                                          exact_num_lineages;
        double                                                          max_time;
        bool                                                            allow_rate_shifts_on_extinct_lineages;
        bool                                                            prune_extinct_lineages;
        bool                                                            condition_on_tip_states;
        bool                                                            condition_on_num_tips;
        bool                                                            condition_on_tree;
        double                                                          NUM_TIME_SLICES;
        std::int64_t                                                    age_check_precision;

    };
    
}

#endif
