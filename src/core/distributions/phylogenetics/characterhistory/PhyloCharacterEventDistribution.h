#ifndef PhyloCharacterEventDistribution_H
#define PhyloCharacterEventDistribution_H

#include "CharacterHistoryDiscrete.h"
#include "MemberObject.h"
//#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class PhyloCharacterEventDistribution : public TypedDistribution< CharacterHistoryDiscrete >, public MemberObject< RbVector<long> >, public MemberObject< RbVector<double> > {
        
    public:
        PhyloCharacterEventDistribution(const RbVector<const TypedDagNode<double>* >& root,
                                        const std::vector<TypedDistribution<double>* >& bd,
                                        const TypedDagNode<Tree>* t,
                                        const TypedDagNode<double> *s,
                                        const std::vector< std::string >& n);                                                                                                           //!< Constructor
        
        virtual                                            ~PhyloCharacterEventDistribution(void);                                                              //!< Virtual destructor
        
        // public member functions
        PhyloCharacterEventDistribution*                    clone(void) const;                                                                                              //!< Create an independent clone
        double                                              computeLnProbability(void);                                                                                     //!< Compute ln prob of current value
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<long> &rv) const;         //!< Map the member methods to internal function calls
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;       //!< Map the member methods to internal function calls
        CharacterHistoryDiscrete&                           getCharacterHistory(void);                                                                                      //!< Get the character histories
        const CharacterHistoryDiscrete&                     getCharacterHistory(void) const;                                                                                //!< Get the character histories
        TypedDistribution<double>*                          getValueDistibution(const std::string& n) const;
        double                                              getRootValue(const std::string& n) const;
        double                                              getRootSpeciationRate(void) const;
        void                                                redrawValue(void);                                                                                              //!< Draw a new random value from distribution
        void                                                setValue(CharacterHistoryDiscrete *v, bool force);
        
    protected:
        // Parameter management functions
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                              //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        double                                              computeNodeProbability(const TopologyNode &n, size_t nIdx);
        
        
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                                //!< Swap a parameter
        
    private:
        
        // helper functions
        void                                                buildRandomBinaryHistory(std::vector<TopologyNode *> &tips);
        double                                              computeStateValue(size_t i, size_t j, double time) const;
        double                                              computeStartValue(size_t i, size_t j) const;
        void                                                initializeBranchHistories(const TopologyNode &n, size_t nIdx);
        
        // members
        RbVector<const TypedDagNode<double>* >              root_values;
        std::vector<TypedDistribution<double>* >            base_distribution;
        const TypedDagNode<Tree>*                           tree;
        const TypedDagNode<double>*                         shift_rate;
        std::vector< std::string >                          names;
        
        size_t                                              num_values_per_event;

        RbVector< RbVector<double> >                        event_values;
//        CharacterHistoryDiscrete                            branch_histories;
                
//        double                                              log_tree_topology_prob;                                                                                         //!< Log-transformed tree topology probability (combinatorial constant).
        
        // only for testing
        bool                                                event_prior_only;
        
    };
    
}

#endif

