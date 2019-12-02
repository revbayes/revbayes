#ifndef ExponentialDemographicFunction_H
#define ExponentialDemographicFunction_H

#include "DemographicFunction.h"
#include "TypedDagNode.h"

#include <iostream>
#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief ExponentialDemographicFunction
     *
     *
     */
    class ExponentialDemographicFunction : public DemographicFunction {
        
    public:
        ExponentialDemographicFunction(const TypedDagNode<double>* N0, const TypedDagNode<double>* N1, const TypedDagNode<double>* t0, const TypedDagNode<double>* t1);                                                                                          //!< Default constructor
        ExponentialDemographicFunction(const ExponentialDemographicFunction &f);                                                                  //!< Copy constructor
        
        // destructor
        virtual                                        ~ExponentialDemographicFunction(void);
        ExponentialDemographicFunction&                 operator=(const ExponentialDemographicFunction &f);                    //!< Assignment operator
        
        // public methods
        
        // pure virtual public methods
        ExponentialDemographicFunction*                 clone(void) const;                                      //!< Clone the ExponentialDemographicFunction
        double                                          getDemographic(double t) const;                         //!< demographic function N(t) at time t
        double                                          getIntegral(double start, double finish) const;         //!< Calculates the integral 1/N(x) dx between start and finish.
        
    protected:
        virtual void                                    swapNodeInternal(const DagNode *oldN, const DagNode *newN);

        
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
