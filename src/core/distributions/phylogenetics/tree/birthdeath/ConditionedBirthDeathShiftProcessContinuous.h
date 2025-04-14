#ifndef ConditionedBirthDeathShiftProcessContinuous_H
#define ConditionedBirthDeathShiftProcessContinuous_H

#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "CharacterHistoryContinuous.h"
#include "MemberObject.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class ConditionedBirthDeathShiftProcessContinuous : public AbstractCharacterHistoryBirthDeathProcess, public MemberObject< RbVector<std::int64_t> >, public MemberObject< RbVector<double> > {
        
    public:
        ConditionedBirthDeathShiftProcessContinuous(const TypedDagNode<double> *a,
                                    const TypedDagNode<double> *root_speciation,
                                    const TypedDagNode<double> *root_extinction,
                                    TypedDistribution<double> *s,
                                    TypedDistribution<double> *e,
                                    const TypedDagNode<double> *ev,
                                    const TypedDagNode<double> *r,
                                    const std::string &cdt,
                                    const std::vector<Taxon> &n);                                                                                                           //!< Constructor
        
        virtual                                            ~ConditionedBirthDeathShiftProcessContinuous(void);                                                              //!< Virtual destructor
        
        // public member functions
        ConditionedBirthDeathShiftProcessContinuous*        clone(void) const;                                                                                              //!< Create an independent clone
        double                                              computeLnProbability(void);                                                                                     //!< Compute ln prob of current value
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<std::int64_t> &rv) const;         //!< Map the member methods to internal function calls
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;       //!< Map the member methods to internal function calls
        CharacterHistoryContinuous&                         getCharacterHistory(void);                                                                                      //!< Get the character histories
        const CharacterHistoryContinuous&                   getCharacterHistory(void) const;                                                                                //!< Get the character histories
        TypedDistribution<double>*                          getExtinctionRateDistibution(void) const;
        TypedDistribution<double>*                          getSpeciationRateDistibution(void) const;
        double                                              getRootExtinctionRate(void) const;
        double                                              getRootSpeciationRate(void) const;
        bool                                                isExtinctionRateConstant(void) const;
        bool                                                isSpeciationRateConstant(void) const;
        void                                                redrawValue(void);                                                                                              //!< Draw a new random value from distribution
        void                                                setValue(Tree *v, bool force);
        
    protected:
        // Parameter management functions
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                              //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        double                                              computeNodeProbability(const TopologyNode &n, size_t nIdx);
        double                                              computeRootLikelihood(void);
        
        
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                                //!< Swap a parameter
        
    private:
        
        // helper functions
        void                                                attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times, double T);
        void                                                buildRandomBinaryHistory(std::vector<TopologyNode *> &tips);
        double                                              computeStateValue(size_t i, size_t j, double time) const;
        double                                              computeStartValue(size_t i, size_t j) const;
        void                                                initializeBranchHistories(const TopologyNode &n, size_t nIdx, size_t idx, double val);
        void                                                simulateTree(void);
        void                                                runMCMC(void);
        double                                              computeBranchProbability(double begin, double end, double sp_rate, double ex_rate, double sh_rate);
        
        // members
        const TypedDagNode<double>*                         root_age;
        const TypedDagNode<double>*                         root_speciation;
        const TypedDagNode<double>*                         root_extinction;
        TypedDistribution<double>*                          base_distribution_speciation;
        TypedDistribution<double>*                          base_distribution_extinction;
        const TypedDagNode<double>*                         shift_rate;
        const TypedDagNode<double>*                         rho;                                                                                                            //!< Sampling probability of each species.
        
        CharacterHistoryContinuous                          branch_histories;
        
        std::string                                         condition;                                                                                                      //!< The condition of the process (none/survival/#taxa).
        size_t                                              num_taxa;
        std::vector<Taxon>                                  taxa;
        
        double                                              log_tree_topology_prob;                                                                                         //!< Log-transformed tree topology probability (combinatorial constant).
        
        
        // only for testing
        bool                                                event_prior_only;
        
    };
    
}

#endif

