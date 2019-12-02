#ifndef LinearDemographicFunction_H
#define LinearDemographicFunction_H

#include "DemographicFunction.h"
#include "TypedDagNode.h"

#include <iostream>
#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief LinearDemographicFunction
     *
     *
     */
    class LinearDemographicFunction : public DemographicFunction {
        
    public:
        LinearDemographicFunction(const TypedDagNode<double>* N0, const TypedDagNode<double>* N1, const TypedDagNode<double>* t0, const TypedDagNode<double>* t1);                                                                                          //!< Default constructor
        LinearDemographicFunction(const LinearDemographicFunction &f);                                                                  //!< Copy constructor
        
        // destructor
        virtual                                        ~LinearDemographicFunction(void);
        LinearDemographicFunction&                      operator=(const LinearDemographicFunction &f);                    //!< Assignment operator
        
        // public methods
        
        // pure virtual public methods
        LinearDemographicFunction*                      clone(void) const;                                      //!< Clone the LinearDemographicFunction
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
    std::ostream&                                       operator<<(std::ostream& o, const LinearDemographicFunction& x);  //!< Overloaded output stream
    
}

#endif /* LinearDemographicFunction_H */
