#ifndef LinearDemographicFunction_H
#define LinearDemographicFunction_H

#include <iostream>

#include "DemographicFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Linear Demographic Function
     *
     * The demographic function related to a linear population growth.
     * It contains four variables:
     *    - theta_recent, the population size at the beginning of the linear period (towards the present)
     *    - theta_ancient, the population size at the end of the linear period (towards the past)
     *    - time_recent, the time of the beginning of the linear period (towards the present)
     *    - time_ancient, the time of the end of the linear period (towards the past)
     */
    class LinearDemographicFunction : public DemographicFunction {
        
    public:
        LinearDemographicFunction(const TypedDagNode<double>* N0, const TypedDagNode<double>* N1, const TypedDagNode<double>* t0, const TypedDagNode<double>* t1);                                                                                          //!< Default constructor
        LinearDemographicFunction(const LinearDemographicFunction &f);                                                                  //!< Copy constructor
        
        // destructor
        virtual                                        ~LinearDemographicFunction(void);                            //!< Destructor
        LinearDemographicFunction&                      operator=(const LinearDemographicFunction &f);              //!< Assignment operator
        
        // public methods
        
        // pure virtual public methods
        LinearDemographicFunction*                      clone(void) const;                                          //!< Clone the LinearDemographicFunction
        double                                          getDemographic(double t) const;                             //!< Returns the demographic function N(t) at time t.
        double                                          getIntegral(double start, double finish) const;             //!< Calculates the integral 1/N(x) dx between start and finish.
        double                                          getWaitingTime(double time, double lambda) const;           //!< Calculates the waiting time until the next coalescent event.
        
    protected:
        virtual void                                    swapNodeInternal(const DagNode *oldN, const DagNode *newN); //!< Internally replacing a DAG node containing a variable with a different one
        
        
    private:
        
        const TypedDagNode<double>*                     theta_ancient;
        const TypedDagNode<double>*                     theta_recent;
        const TypedDagNode<double>*                     time_ancient;
        const TypedDagNode<double>*                     time_recent;
        
        
    };
    
    // Global functions using the class
    std::ostream&                                       operator<<(std::ostream& o, const LinearDemographicFunction& x);  //!< Overloaded output stream
    
}

#endif /* LinearDemographicFunction_H */
