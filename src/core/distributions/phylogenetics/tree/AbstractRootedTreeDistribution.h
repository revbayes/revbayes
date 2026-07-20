#ifndef AbstractRootedTreeDistribution_H
#define AbstractRootedTreeDistribution_H

#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * Constant rate Birth-Death process.
     *
     */
    class AbstractRootedTreeDistribution : public TypedDistribution<Tree> {
        
    public:
        AbstractRootedTreeDistribution(const TypedDagNode<double> *ra, const std::vector<Taxon> &tn, bool uo, Tree *t );
        AbstractRootedTreeDistribution(const AbstractRootedTreeDistribution& d );

        
        virtual ~AbstractRootedTreeDistribution(void);
        
        // overloaded operators
        AbstractRootedTreeDistribution&                     operator=(const AbstractRootedTreeDistribution& d);                                                            //!< Assignment

        // pure virtual member functions
        virtual AbstractRootedTreeDistribution*             clone(void) const = 0;                                                                              //!< Create an independent clone
        
        
        // public member functions you may want to override
        double                                              computeLnProbability(void);                                                                         //!< Compute the log-transformed probability of the current value.
        virtual void                                        redrawValue(SimulationCondition c);                                                                 //!< Draw a new random value from the distribution
        virtual void                                        redrawValue(void);                                                                                  //!< Draw a new random value from the distribution

        virtual void                                        setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        size_t                                              getNumberOfTaxa(void) const;
        double                                              getOriginAge(void) const;
        virtual double                                      getRootAge(void) const;
        const std::vector<Taxon>&                           getTaxa(void) const;
        virtual void                                        simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn);
        virtual double                                      simulateCladeAge(size_t n, double origin, double present, double min, bool alwaysReturn) const;

    protected:
        // pure virtual helper functions
        virtual double                                      computeLnProbabilityDivergenceTimes(void) const = 0;                                                //!< Compute the log-transformed probability of the current value.

        virtual bool                                        isLnProbabilityNonZero(void);
        virtual double                                      simulateDivergenceTime(double origin, double present) const = 0;                                    //!< Simulate n speciation events.
        virtual std::vector<double>                         simulateDivergenceTimes(size_t n, double origin, double end, double present, bool alwaysReturn) const = 0;             //!< Simulate n speciation events.
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                  //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // Parameter management functions. You need to override both if you have additional parameters
        virtual void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        
        
        // helper functions
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        virtual double                                      lnProbTreeShape(void) const;
        void                                                recomputeDivergenceTimesSinceOrigin(void) const;                                                    //!< Extract the divergence times from the tree.
        int                                                 diversity(double t);                                                                                //!< Diversity at time t.
        std::vector<double>                                 getAgesOfInternalNodesFromMostRecentSample(void) const;                                             //!< Get the ages of all internal nodes since the time of the most recent tip age.
        std::vector<double>                                 getAgesOfTipsFromMostRecentSample(void) const;                                                      //!< Get the ages of all tip nodes since the time of the most recent tip age.

        double                                              simulateNextAge(size_t n, double origin, double present, double min, bool alwaysReturn) const;
        void                                                simulateTree(bool alwaysReturn);
        
        // members
        mutable std::vector<double>                         divergence_times;                                                                                   //!< Taxon names that will be attached to new simulated trees.

        const TypedDagNode<double>*                         process_age;                                                                                        //!< Time since the start of the process.
        std::vector<Taxon>                                  taxa;                                                                                               //!< Taxon names that will be attached to new simulated trees.
        bool                                                use_origin;
        Tree*                                               starting_tree;
    };
    
}

#endif
