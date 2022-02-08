#ifndef PiecewiseConstantHeterochronousCoalescent_H
#define PiecewiseConstantHeterochronousCoalescent_H

#include "AbstractCoalescent.h"
#include "RbVector.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class Clade;
    
    /**
     * @brief Piecewise-Constant Heterochronous Coalescent Process.
     *
     * This process is an extension to the piecewise constant coalescent to include serially sampled lineages.
     * This process has intervals where population sizes are constant within each interval.
     *
     * The process can have one or two parameters:
     * @param NE    The population sizes for each interval
     * @param internvalStarts tTe start of a new interval (0 is implicitely assumed). The start of a new interval denotes the end of the previous interval.
     *
     */
    class PiecewiseConstantHeterochronousCoalescent : public AbstractCoalescent {
        
    public:
        
        enum METHOD_TYPES { EVENTS, SPECIFIED, UNIFORM };
        
        PiecewiseConstantHeterochronousCoalescent(const TypedDagNode<RbVector<double> > *N, const TypedDagNode<RbVector<double> > *i, const std::vector<Taxon> &tn, const std::vector<Clade> &c);
        virtual                                            ~PiecewiseConstantHeterochronousCoalescent(void);                                                                    //!< Virtual destructor
        
        // public member functions
        PiecewiseConstantHeterochronousCoalescent*          clone(void) const;                                                                                  //!< Create an independent clone
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
        // derived helper functions
        double                                              computeLnProbabilityTimes(void) const;                                                          //!< Compute the log-transformed probability of the current value.
        std::vector<double>                                 simulateCoalescentAges(size_t n) const;                                                         //!< Simulate n coalescent events.
        
        
    private:
        
        
        // members
        const TypedDagNode<RbVector<double> >*              Nes;                            //!< The population size for each interval
        const TypedDagNode<RbVector<double> >*              intervalStarts;                 //!< The time of a new interval start time.
        
    };
    
}

#endif
