#ifndef PiecewiseCoalescent_H
#define PiecewiseCoalescent_H

#include "AbstractCoalescent.h"
#include "DemographicFunction.h"
#include "ConstantDemographicFunction.h"
#include "LinearDemographicFunction.h"
#include "RbVector.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class Clade;
    
    /**
     * @brief Piecewise-constant population size coalescent process.
     *
     *
     * The piecewise-constant population size coalescent process is an extension to the constant
     * population size coalescent process. Instead of having a constant population size over the duration of the process, the process is split into intervals where the population size within each interval is constant.
     *  The process can have one or two parameters:
     * @param NEs the population sizes
     * @param internvalStarts the start time of a new interval (0 is implicitely assumed for the first interval)
     *
     *
     *
     *
     */
    class PiecewiseCoalescent : public AbstractCoalescent, public MemberObject< RbVector<double> > {
        
    public:
        
        enum METHOD_TYPES { EVENTS, SPECIFIED };
        enum DEMOGRAPHY_FUNCTION_TYPES { CONSTANT, LINEAR };

        
        PiecewiseCoalescent(const TypedDagNode<RbVector<double> > *N, const TypedDagNode<RbVector<double> > *i, const TypedDagNode<RbVector<std::int64_t> > *n_events_pi, METHOD_TYPES meth, DEMOGRAPHY_FUNCTION_TYPES dem, const std::vector<Taxon> &tn, const std::vector<Clade> &c);
        virtual                                            ~PiecewiseCoalescent(void);                                                                    //!< Virtual destructor
        
        // public member functions
        PiecewiseCoalescent*                                clone(void) const;                                                                                  //!< Create an independent clone
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls

    protected:
        // Parameter management functions
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                //!< Swap a parameter
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);

        // derived helper functions
        double                                              computeLnProbabilityTimes(void) const;                                                          //!< Compute the log-transformed probability of the current value.
        std::vector<double>                                 simulateCoalescentAges(size_t n) const;
        
        
    private:
        
        void                                                updateIntervals(void) const;
        double                                              getDemographic(double event_age, double index) const;
        double                                              getIntegral(double last_age, double event_age, double index) const;
        double                                              getWaitingTime(double age, double rv, size_t index) const;

        // members
        const TypedDagNode<RbVector<double> >*              Nes;                                    //!< A pointer for the population sizes for each interval
        const TypedDagNode<RbVector<double> >*              interval_change_points_var;
        const TypedDagNode<RbVector<std::int64_t> >*                number_events_per_interval;
        mutable RbVector<double>                            interval_change_points;
        mutable RbVector<double>                            pop_sizes;                              //!< The population sizes for each interval
        METHOD_TYPES                                        interval_method;                        //!< The method of specifying coalescent intervals
        DEMOGRAPHY_FUNCTION_TYPES                           demographic_function_var;

    };
    
}

#endif
