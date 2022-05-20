#ifndef DirichletTimeTreeDistribution_H
#define DirichletTimeTreeDistribution_H

#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class DirichletTimeTreeDistribution : public TypedDistribution<Tree> {
        
    public:
        DirichletTimeTreeDistribution(const TypedDagNode<double> *r, const TypedDagNode< RbVector<double> > *a, const std::vector<Taxon> &n);                                                                                  //!< Constructor

        virtual                                            ~DirichletTimeTreeDistribution(void);                          //!< Virtual destructor
        
        // public member functions
        DirichletTimeTreeDistribution*                      clone(void) const;                                          //!< Create an independent clone
        double                                              computeLnProbability(void);                                 //!< Compute ln prob of current value
        void                                                redrawValue(void);                                          //!< Draw a new random value from distribution

        
    protected:
        // Parameter management functions
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                      //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);

        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // helper functions
        void                                                attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times, double T);
        void                                                buildRandomBinaryHistory(std::vector<TopologyNode *> &tips);
        void                                                simulateTree(void);
        
        // members
        const TypedDagNode<double>*                         root_age;
        const TypedDagNode< RbVector<double> >*             alpha;
        size_t                                              num_taxa;
        std::vector<Taxon>                                  taxa;
    };
    
}

#endif
