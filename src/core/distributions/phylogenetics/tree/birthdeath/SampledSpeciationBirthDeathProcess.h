#ifndef SampledSpeciationBirthDeathProcess_h
#define SampledSpeciationBirthDeathProcess_h

#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "CharacterHistoryDiscrete.h"
#include "MemberObject.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class SampledSpeciationBirthDeathProcess : public AbstractCharacterHistoryBirthDeathProcess, public MemberObject< RbVector<long> >, public MemberObject<RbVector<RbVector<double> > > {
        
    public:
        SampledSpeciationBirthDeathProcess(const TypedDagNode<double> *a, const TypedDagNode<double> *s, const TypedDagNode<double> *e, const TypedDagNode<double> *r, const std::vector<Taxon> &n);                                                                                  //!< Constructor
        
        virtual                                             ~SampledSpeciationBirthDeathProcess(void);                          //!< Virtual destructor
        
        // public member functions
        SampledSpeciationBirthDeathProcess*                 clone(void) const;                                          //!< Create an independent clone
        double                                              computeLnProbability(void);                                 //!< Compute ln prob of current value
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<long> &rv) const;     //!< Map the member methods to internal function calls
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<double> > &rv) const;     //!< Map the member methods to internal function calls

        CharacterHistoryDiscrete&                           getCharacterHistory(void);                                  //!< Get the character histories
        const CharacterHistoryDiscrete&                     getCharacterHistory(void) const;                                  //!< Get the character histories
        void                                                redrawValue(void);                                          //!< Draw a new random value from distribution
        void                                                setValue(Tree *v, bool force);
        
    protected:
        // Parameter management functions
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter);                                      //!< get affected nodes
        virtual void                                        keepSpecialization(DagNode* affecter);
        virtual void                                        restoreSpecialization(DagNode *restorer);
        virtual void                                        touchSpecialization(DagNode *toucher, bool touchAll);
        
        void                                                computeNodeProbability(const TopologyNode &n, size_t nIdx);
        double                                              computeRootLikelihood(void);
        
        
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // helper functions
        void                                                assignNodes(TopologyNode* n, size_t& tip_index, size_t& int_index);
        void                                                assertParentChildEdges(TopologyNode* n);
        double                                              computeLineageUnsampledByPresentProbability(double t, double T);
//        double                                              computeLnProbabilityTaxa(void);
//        double                                              recursivelyComputeLnProbabilityTaxa(TopologyNode* nd);
        void                                                getLineagesAtAge(TopologyNode* n, std::vector<TopologyNode*>& nodes, double t);
        void                                                simulateEventsForTreeAdHoc(void);
        size_t                                              simulateEvent(TopologyNode* n, std::set<TopologyNode*>& nodes, std::vector<double>& unsampledLineages, double time, double maxTime);
        void                                                simulateTree(void);
        void                                                simulateUnsampledLineages(Tree* t, std::vector<double> ages);
        void                                                simulateUnsampledRootLineages(Tree* t, double rootAge);


        
        // members
        const TypedDagNode<double>*                         root_age;
        const TypedDagNode<double>*                         speciation;
        const TypedDagNode<double>*                         extinction;
        const TypedDagNode<double>*                         rho;                                                //!< Sampling probability of each species.
        
        CharacterHistoryDiscrete                            branch_histories;
        
        size_t                                              num_taxa;
        std::vector<Taxon>                                  taxa;
        
        
        std::vector<size_t>                                 activeLikelihood;
        std::vector<std::vector<double> >                   storedLikelihood;
        mutable std::vector<bool>                           changed_nodes;
        mutable std::vector<bool>                           dirty_nodes;
        mutable std::vector<std::vector<double> >           scalingFactors;
        mutable double                                      totalScaling;
        
        double                                              logTreeTopologyProb;                                                                                //!< Log-transformed tree topology probability (combinatorial constant).
        
    };
    
}

#endif /* SampledSpeciationBirthDeathProcess_h */
