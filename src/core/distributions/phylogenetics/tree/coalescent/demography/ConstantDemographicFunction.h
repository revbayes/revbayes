#ifndef ConstantDemographicFunction_H
#define ConstantDemographicFunction_H

#include <iostream>

#include "DemographicFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Constant demographic function.
     *
     * The constant demographic function is the simplest available demographic function.
     * It contains only a single variable, theta, the population size.
     */
    class ConstantDemographicFunction : public DemographicFunction {
        
    public:
        ConstantDemographicFunction(const TypedDagNode<double>* t);                                             //!< Default constructor
        ConstantDemographicFunction(const ConstantDemographicFunction &f);                                      //!< Copy constructor
        
        // destructor
        virtual                                        ~ConstantDemographicFunction(void);                      //!< Destructor

        ConstantDemographicFunction&                    operator=(const ConstantDemographicFunction &f);        //!< Assignment operator

        // public methods
        
        // pure virtual public methods
        ConstantDemographicFunction*                    clone(void) const;                                      //!< Clone the ConstantDemographicFunction
        double                                          getDemographic(double t) const;                         //!< Returns the demographic function N(t) at time t.
        double                                          getIntegral(double start, double finish) const;         //!< Calculates the integral 1/N(x) dx between start and finish.
        double                                          getWaitingTime(double time, double lambda, double ploidy=1.0) const;       //!< Calculates the waiting time until the next coalescent event.
        
    protected:
        virtual void                                    swapNodeInternal(const DagNode *oldN, const DagNode *newN); //!< Internally replacing a DAG node containing a variable with a different one

        
    private:
        
        const TypedDagNode<double>*                     theta;

        
    };
    
    // Global functions using the class
    std::ostream&                                       operator<<(std::ostream& o, const ConstantDemographicFunction& x);  //!< Overloaded output stream
    
}

#endif /* ConstantDemographicFunction_H */
