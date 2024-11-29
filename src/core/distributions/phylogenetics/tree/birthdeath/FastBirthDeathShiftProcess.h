#ifndef FastBirthDeathShiftProcess_H
#define FastBirthDeathShiftProcess_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TreeDiscreteCharacterData.h"
#include "RateMatrix.h"
#include "RateMatrix_JC.h"
#include "Simplex.h"
#include "BDS_ODE.h"
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
    class FastBirthDeathShiftProcess : public TypedDistribution<Tree>, public TreeChangeEventListener, public MemberObject< RbVector<long> >, public MemberObject< RbVector<double> > {
        
    public:
        FastBirthDeathShiftProcess(const TypedDagNode<double> *root,
                                  const TypedDagNode<RbVector<double> >* s,
                                  const TypedDagNode<RbVector<double> >* m,
                                  const TypedDagNode<double> * rsp,
                                  const TypedDagNode<double> * rext,
                                  const TypedDagNode<RateGenerator>* q,
                                  const TypedDagNode<double>* r,
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
                                  bool allow_shifts_extinct);
        
        // pure virtual member functions
        virtual FastBirthDeathShiftProcess*                             clone(void) const;
        virtual                                                         ~FastBirthDeathShiftProcess(void);                                                              //!< Virtual destructor

        double                                                          computeLnProbability(void);
        void                                                            fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        const AbstractHomologousDiscreteCharacterData&                  getCharacterData() const;
        double                                                          getOriginAge(void) const;
        std::vector<double>                                             getAverageExtinctionRatePerBranch(void) const;
        std::vector<double>                                             getAverageSpeciationRatePerBranch(void) const;
        std::vector<long>                                               getNumberOfShiftEventsPerBranch(void) const;
        std::vector<double>                                             getTimeInStates(void) const;
        double                                                          getRootAge(void) const;
        virtual void                                                    redrawValue(void);
        void                                                            setSampleCharacterHistory(bool sample_history);                                                     //!< Set whether or not we are sampling the character history along branches.
        void                                                            setSamplingFraction(const TypedDagNode< double > *r);
        void                                                            setSamplingFraction(const TypedDagNode< RbVector<double> > *r);
        void                                                            setSerialSamplingRates(const TypedDagNode< RbVector<double> > *r);
        virtual void                                                    setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
        void                                                            drawJointConditionalAncestralStates(std::vector<size_t>& startStates, std::vector<size_t>& endStates);
        void                                                            recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<size_t>& startStates, std::vector<size_t>& endStates);
        void                                                            numericallyIntegrateProcess(std::vector< double > &likelihoods, double begin_age, double end_age, bool use_backward, bool extinction_only) const; //!< Wrapper function for the ODE time stepper function.
        void                                                            resizeVectors(size_t num_nodes);
        
    protected:
        
        double                                                          getEventRate(void) const;
        const RateGenerator&                                            getEventRateMatrix(void) const;
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
        void                                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<long> &rv) const;     //!< Map the member methods to internal function calls
        RevLanguage::RevPtr<RevLanguage::RevVariable>                   executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found);
        
        // helper functions
        void                                                            buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        std::vector<double>                                             pExtinction(double start, double end) const;                                                        //!< Compute the probability of extinction of the process (without incomplete taxon sampling).
        virtual double                                                  pSurvival(double start, double end) const;                                                          //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double                                                          pSurvival(double start, double end, bool speciation) const;                                                          //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        void                                                            recursivelyFlagNodeDirty(const TopologyNode& n);
        bool                                                            simulateTree(size_t attempts = 0);
        bool                                                            simulateTreeConditionedOnTips(size_t attempts = 0);
        std::vector<double>                                             calculateTotalAnageneticRatePerState(void) const;
        std::vector<double>                                             calculateTotalSpeciationRatePerState(void) const;
        void                                                            computeNodeProbability(const TopologyNode &n, size_t nIdx) const;
        double                                                          computeRootLikelihood() const;
        
        // members
        std::string                                                     condition;                                                                                          //!< The condition of the process (none/survival/#taxa).
        double                                                          dt;                                                                                                 //!< The size of the time slices used by the ODE for numerical integration.
        std::vector<bool>                                               active_likelihood;
        mutable std::vector<bool>                                       changed_nodes;
        mutable std::vector<bool>                                       dirty_nodes;
        mutable std::vector<std::vector<std::vector<double> > >         node_partial_likelihoods;
        mutable std::map<size_t, std::vector<std::vector<double> > >    branch_partial_likelihoods;
        mutable std::vector<std::vector<double> >                       extinction_probabilities;
        size_t                                                          num_states;
        mutable std::vector<std::vector<double> >                       scaling_factors;
        bool                                                            use_origin;
        bool                                                            sample_character_history;                                                                           //!< are we sampling the character history along branches?
        std::vector<double>                                             average_speciation;
        std::vector<double>                                             average_extinction;
        std::vector<long>                                               num_shift_events;
        std::vector<double>                                             time_in_states;
        std::string                                                     simmap;
        
        // parameters
        const TypedDagNode<double>*                                     process_age;                                                                                           //!< Time since the origin.
        const TypedDagNode<RbVector<double> >*                          mu;
        const TypedDagNode<RbVector<double> >*                          lambda;
        const TypedDagNode<double>*                                     alpha;
        const TypedDagNode<double>*                                     beta;
        const TypedDagNode<Simplex >*                                   pi;                                                                                                 //!< The root frequencies (probabilities of the root states).
        const TypedDagNode<RateGenerator>*                              Q;
        const TypedDagNode<double>*                                     rate;                                                                                               //!< Sampling probability of each species.
        const TypedDagNode<double>*                                     rho;                                                                                                //!< Sampling probability of each species.
        const TypedDagNode<RbVector<double> >*                          rho_per_state;                                                                                                //!< Sampling probability of each species.

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
    };
    
}

#endif
