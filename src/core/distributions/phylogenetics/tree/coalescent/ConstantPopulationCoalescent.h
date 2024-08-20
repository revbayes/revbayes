#ifndef ConstantPopulationCoalescent_H
#define ConstantPopulationCoalescent_H

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
     * It contains only a single additional parameter from the abstract class.
     *
     * @param Ne the population size
     *
     * @copydoc RevBayesCore::AbstractCoalescent
     *
     * @see RevBayesCore::AbstractCoalescent for the parent class
     *
     *
     */
    class ConstantPopulationCoalescent : public AbstractCoalescent {
        
    public:
        ConstantPopulationCoalescent(const TypedDagNode<double> *N, const std::vector<Taxon> &tn, const std::vector<Clade> &c);
        virtual                                            ~ConstantPopulationCoalescent(void);                                                             //!< Virtual destructor
        
        // public member functions
        ConstantPopulationCoalescent*                       clone(void) const;                                                                              //!< Create an independent clone

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                //!< Swap a parameter
        
        // derived helper functions
        double                                              computeLnProbabilityTimes(void) const;                                                          //!< Compute the log-transformed probability of the current value.
        std::vector<double>                                 simulateCoalescentAges(size_t n) const;                                                         //!< Simulate n coalescent events.

        
    private:
        
        
        // members
        const TypedDagNode<double>*                         Ne; //!< The effective population size
        
    };
    
}

#endif
