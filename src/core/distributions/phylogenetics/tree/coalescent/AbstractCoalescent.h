#ifndef AbstractCoalescent_H
#define AbstractCoalescent_H

#include <cstddef>
#include <vector>

#include "Taxon.h"
#include "Tree.h"
#include "TypedDistribution.h"
#include "Clade.h"

namespace RevBayesCore {
    
    class TopologyNode;
    
    /**
     * @brief Declaration of the abstract coalescent process class.
     *
     *This file contains the declaration of the random variable class for any coalescent process.
     */
    class AbstractCoalescent : public TypedDistribution<Tree> {
        
    public:
        AbstractCoalescent(const std::vector<Taxon> &tn, const std::vector<Clade> &c);
        
        // pure virtual member functions
        virtual AbstractCoalescent*                         clone(void) const = 0;                                                                              //!< Create an independent clone
        
        
        // public member functions you may want to override
        double                                              computeLnProbability(void);                                                                         //!< Compute the log-transformed probability of the current value.
        virtual void                                        redrawValue(void);                                                                                  //!< Draw a new random value from the distribution
        
        
    protected:
        
        // pure virtual helper functions
        virtual double                                      computeLnProbabilityTimes(void) const = 0;                                                          //!< Compute the log-transformed probability of the current value.
        virtual std::vector<double>                         simulateCoalescentAges(size_t n) const = 0;                                                         //!< Simulate n coalescent events.
        
        // helper functions
        void                                                attachAges(Tree *psi, std::vector<TopologyNode *> &tips, size_t index,
                                                                        const std::vector<double> &a);
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        void                                                buildHeterochronousRandomBinaryTree(Tree *psi, std::vector<TopologyNode*> &active, const std::vector<double> &ages);
        bool                                                matchesConstraints(void);
        void                                                simulateTree(void);
        void                                                simulateHeterochronousTree(void);                                                                   //!< Simulates a heterochronus coalescent tree.
        
        // members
        std::vector<Clade>                                  constraints;                                                                                        //!< Topological constrains.
        size_t                                              num_taxa;                                                                                           //!< Number of taxa (needed for correct initialization).
        std::vector<Taxon>                                  taxa;                                                                                               //!< Taxon names that will be attached to new simulated trees.
        double                                              logTreeTopologyProb;                                                                                //!< Log-transformed tree topology probability (combinatorial constant).
        
    };
    
}

#endif
