#ifndef HeterochronousCoalescent_H
#define HeterochronousCoalescent_H

#include "AbstractCoalescent.h"
#include "DemographicFunction.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class Clade;
    
    
    
    /**
     * @brief Heterochronous Coalescent Process.
     *
     * This process involves coalescence of serially sampled lineages. The process is split into intervals where the demography dynamics can vary between intervals.
     */
    class HeterochronousCoalescent : public AbstractCoalescent {
        
    public:
        HeterochronousCoalescent(const TypedDagNode< RbVector<double> > *iv, const RbVector< DemographicFunction > &df, const std::vector<Taxon> &tn, const std::vector<Clade> &c);
//        HeterochronousCoalescent(const HeterochronousCoalescent &d);
        virtual                                            ~HeterochronousCoalescent(void);                                                                 //!< Virtual destructor
        
//        HeterochronousCoalescent&                           operator=(const HeterochronousCoalescent &d);

        // public member functions
        HeterochronousCoalescent*                           clone(void) const;                                                                              //!< Create an independent clone
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                //!< Swap a parameter
        
        // derived helper functions
        double                                              computeLnProbabilityTimes(void) const;                                                          //!< Compute the log-transformed probability of the current value.
        std::vector<double>                                 simulateCoalescentAges(size_t n) const;                                                         //!< Simulate n coalescent events.
        
        
    private:
        
        enum EVENT_TYPE { COALESCENT, SERIAL_SAMPLE, DEMOGRAPHIC_MODEL_CHANGE };
        
        // members
        const TypedDagNode< RbVector<double> >*             intervals; //!<The start times for intervals
        RbVector< DemographicFunction >                     demographies; //!< a vector of functions that model how the demogrpahy changes over the course of that interval
    };
    
}

#endif
