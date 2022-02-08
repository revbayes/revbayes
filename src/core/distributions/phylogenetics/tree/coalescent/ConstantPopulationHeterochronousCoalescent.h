#ifndef ConstantPopulationHeterochronousCoalescent_H
#define ConstantPopulationHeterochronousCoalescent_H

#include "AbstractCoalescent.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class Clade;
    
    
    
    /**
     * @brief Constant population size coalescent process.
     *
     *
     * The constant population size coalescent process is the simplest available coalescent process.
     * It contains only a single parameter NE, the population size.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-17, version 1.0
     *
     */
    class ConstantPopulationHeterochronousCoalescent : public AbstractCoalescent {
        
    public:
        ConstantPopulationHeterochronousCoalescent(const TypedDagNode<double> *N, const std::vector<Taxon> &tn, const std::vector<Clade> &c);
        virtual                                            ~ConstantPopulationHeterochronousCoalescent(void);                                               //!< Virtual destructor
        
        // public member functions
        ConstantPopulationHeterochronousCoalescent*         clone(void) const;                                                                              //!< Create an independent clone

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                //!< Swap a parameter
        
        // derived helper functions
        double                                              computeLnProbabilityTimes(void) const;                                                          //!< Compute the log-transformed probability of the current value.
        std::vector<double>                                 simulateCoalescentAges(size_t n) const;                                                         //!< Simulate n coalescent events.

        
    private:
        
        
        // members
        const TypedDagNode<double>*                         Ne;
        
    };
    
}

#endif
