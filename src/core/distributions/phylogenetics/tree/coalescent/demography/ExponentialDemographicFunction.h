#ifndef ExponentialDemographicFunction_H
#define ExponentialDemographicFunction_H

#include <iostream>

#include "DemographicFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Exponential demographic function
     *
     * The demographic function related to an exponential population growth.
     * It contains four variables:
     *    - theta_recent, the population size at the beginning of the exponential period (towards the present)
     *    - theta_ancient, the population size at the end of the exponential period (towards the past)
     *    - time_recent, the time of the beginning of the exponential period (towards the present)
     *    - time_ancient, the time of the end of the exponential period (towards the past)
     */
    class ExponentialDemographicFunction : public DemographicFunction {
        
    public:
        ExponentialDemographicFunction(const TypedDagNode<double>* N0, const TypedDagNode<double>* N1, const TypedDagNode<double>* t0, const TypedDagNode<double>* t1);                                                                                          //!< Default constructor
        ExponentialDemographicFunction(const ExponentialDemographicFunction &f);                                    //!< Copy constructor
        
        // destructor
        virtual                                        ~ExponentialDemographicFunction(void);                       //!< Destructor
        ExponentialDemographicFunction&                 operator=(const ExponentialDemographicFunction &f);         //!< Assignment operator
        
        // public methods
        
        // pure virtual public methods
        ExponentialDemographicFunction*                 clone(void) const;                                          //!< Clone the ExponentialDemographicFunction
        double                                          getDemographic(double t) const;                             //!< Returns the demographic function N(t) at time t.
        double                                          getIntegral(double start, double finish) const;             //!< Calculates the integral 1/N(x) dx between start and finish.
        double                                          getWaitingTime(double time, double lambda, double ploidy=1.0) const;           //!< Calculates the waiting time until the next coalescent event.
        
    protected:
        virtual void                                    swapNodeInternal(const DagNode *oldN, const DagNode *newN); //!< Internally replacing a DAG node containing a variable with a different one

        
    private:
        
        const TypedDagNode<double>*                     theta_ancient;
        const TypedDagNode<double>*                     theta_recent;
        const TypedDagNode<double>*                     time_ancient;
        const TypedDagNode<double>*                     time_recent;

        
    };
    
    // Global functions using the class
    std::ostream&                                       operator<<(std::ostream& o, const ExponentialDemographicFunction& x);  //!< Overloaded output stream
    
}

#endif /* ExponentialDemographicFunction_H */
