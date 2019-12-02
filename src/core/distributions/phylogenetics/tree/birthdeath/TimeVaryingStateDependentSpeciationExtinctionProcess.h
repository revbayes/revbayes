#ifndef TimeVaryingStateDependentSpeciationExtinctionProcess_H
#define TimeVaryingStateDependentSpeciationExtinctionProcess_H

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
#include "RbVector.h"
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
    class TimeVaryingStateDependentSpeciationExtinctionProcess : public TypedDistribution<Tree>, public TreeChangeEventListener, public MemberObject< RbVector<double> > {
        
    public:
        TimeVaryingStateDependentSpeciationExtinctionProcess(const TypedDagNode<double> *root,
                                                  const TypedDagNode<RbVector< RbVector<double> > >* m,
                                                  const TypedDagNode<RbVector< RateGenerator > >* q,
                                                  const TypedDagNode<RbVector< double > >* et,
                                                  const TypedDagNode<double>* r,
                                                  const TypedDagNode<Simplex>* p,
                                                  const TypedDagNode<double> *rh,
                                                  const std::string &cdt,
                                                  bool uo,
                                                  size_t min_lineages,
                                                  size_t max_lineages,
                                                  bool prune);
        
        // pure virtual member functions
        virtual TimeVaryingStateDependentSpeciationExtinctionProcess*   clone(void) const;
        virtual                                                         ~TimeVaryingStateDependentSpeciationExtinctionProcess(void);                                                              //!< Virtual destructor
        
        double                                                          computeLnProbability(void);
        void                                                            fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        const AbstractHomologousDiscreteCharacterData&                  getCharacterData() const;
        double                                                          getOriginAge(void) const;
        std::vector<double>                                             getTimeInStates(void) const;
        double                                                          getRootAge(void) const;
        virtual void                                                    redrawValue(void);
        void                                                            setCladogenesisMatrix(const TypedDagNode< CladogeneticSpeciationRateMatrix > *r);
        void                                                            setSerialSamplingRates(const TypedDagNode< RbVector<RbVector<double> > > *r);
        void                                                            setSampleCharacterHistory(bool sample_history);                                                     //!< Set whether or not we are sampling the character history along branches.
        void                                                            setSpeciationRates(const TypedDagNode< RbVector< RbVector<double> > > *r);
        void                                                            setNumberOfTimeSlices(double n);                                                                    //!< Set the number of time slices for the numerical ODE.
        virtual void                                                    setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
        void                                                            drawJointConditionalAncestralStates(std::vector<size_t>& startStates, std::vector<size_t>& endStates);
        void                                                            recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<size_t>& startStates, std::vector<size_t>& endStates);
        void                                                            drawStochasticCharacterMap(std::vector<std::string*>& character_histories);
        void                                                            recursivelyDrawStochasticCharacterMap(const TopologyNode &node, size_t start_state, std::vector<std::string*>& character_histories);
        void                                                            numericallyIntegrateProcess(state_type &likelihoods, double begin_age, double end_age, bool use_backward, bool extinction_only) const; //!< Wrapper function for the ODE time stepper function.
        void                                                            resizeVectors(size_t num_nodes);
        
    protected:
        
        double                                                          getEventRate(void) const;
        const RateGenerator&                                            getEventRateMatrix(double a=-1.0) const;
        std::vector<double>                                             getRootFrequencies(void) const;
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                    getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter);                                  //!< get affected nodes
        virtual void                                                    keepSpecialization(DagNode* affecter);
        virtual void                                                    restoreSpecialization(DagNode *restorer);
        virtual void                                                    touchSpecialization(DagNode *toucher, bool touchAll);
        
        double                                                          lnProbTreeShape(void) const;
        
        // Parameter management functions. You need to override both if you have additional parameters
        virtual void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        void                                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;
        RevLanguage::RevPtr<RevLanguage::RevVariable>                   executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found);
        
        // helper functions
        void                                                            buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        std::vector<double>                                             calculateTotalAnageneticRatePerState(void);
        std::vector<double>                                             calculateTotalSpeciationRatePerState(double a=-1.0);
        std::vector<double>                                             calculateExtinctionRatePerState(double a);
        size_t                                                          computeEpochIndex(double a) const;
        double                                                          computeEpochEnd(size_t i) const;
        void                                                            computeNodeProbability(const TopologyNode &n, size_t nIdx) const;
        double                                                          computeRootLikelihood() const;
        std::vector<double>                                             pExtinction(double start, double end) const;                                                        //!< Compute the probability of extinction of the process (without incomplete taxon sampling).
        virtual double                                                  pSurvival(double start, double end) const;                                                          //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        void                                                            recursivelyFlagNodeDirty(const TopologyNode& n);
        bool                                                            simulateTree(size_t attempts = 0);
        
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
        bool                                                            use_cladogenetic_events;                                                                            //!< do we use the speciation rates from the cladogenetic event map?
        bool                                                            use_origin;
        bool                                                            sample_character_history;                                                                           //!< are we sampling the character history along branches?
        std::vector<double>                                             average_speciation;
        std::vector<double>                                             average_extinction;
        std::vector<double>                                             time_in_states;
        std::string                                                     simmap;
        
        // parameters
        const TypedDagNode< CladogeneticSpeciationRateMatrix >*         cladogenesis_matrix;
        const TypedDagNode<double>*                                     process_age;                                                                                           //!< Time since the origin.
        const TypedDagNode<RbVector<RbVector<double> > >*               mu;
        const TypedDagNode<RbVector<RbVector<double> > >*               lambda;
        const TypedDagNode<RbVector<RbVector<double> > >*               phi;
        const TypedDagNode<RbVector<double> >*                          epoch_times;
        const TypedDagNode<Simplex >*                                   pi;                                                                                                 //!< The root frequencies (probabilities of the root states).
        const TypedDagNode<RbVector<RateGenerator> >*                   Q;
        const TypedDagNode<double>*                                     rate;                                                                                               //!< Sampling probability of each species.
        const TypedDagNode<double>*                                     rho;                                                                                                //!< Sampling probability of each species.
        
        RateMatrix_JC                                                   Q_default;
        size_t                                                          min_num_lineages;
        size_t                                                          max_num_lineages;
        bool                                                            prune_extinct_lineages;
        double                                                          NUM_TIME_SLICES;
    };
    
}

#endif
