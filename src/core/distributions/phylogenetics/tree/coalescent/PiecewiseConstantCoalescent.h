#ifndef PiecewiseConstantCoalescent_H
#define PiecewiseConstantCoalescent_H

#include "AbstractCoalescent.h"
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
    class PiecewiseConstantCoalescent : public AbstractCoalescent, public MemberObject< RbVector<double> > {
        
    public:
        
        enum METHOD_TYPES { EVENTS, SPECIFIED, UNIFORM };

        
        PiecewiseConstantCoalescent(const TypedDagNode<RbVector<double> > *N, const TypedDagNode<RbVector<double> > *i, METHOD_TYPES meth, const std::vector<Taxon> &tn, const std::vector<Clade> &c);
        virtual                                            ~PiecewiseConstantCoalescent(void);                                                                    //!< Virtual destructor
        
        // public member functions
        PiecewiseConstantCoalescent*                        clone(void) const;                                                                                  //!< Create an independent clone
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
        
        // members
        const TypedDagNode<RbVector<double> >*              Nes;                                    //!< A pointer for the population sizes for each interval
        const TypedDagNode<RbVector<double> >*              interval_change_points_var;
        mutable RbVector<double>                            interval_change_points;
        mutable RbVector<double>                            pop_sizes;                              //!< The population sizes for each interval
        METHOD_TYPES                                        interval_method;                        //!< The method of specifying coalescent intervals

    };
    
}

#endif
