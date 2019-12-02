#ifndef ConstantDemographicFunction_H
#define ConstantDemographicFunction_H

#include <iostream>

#include "DemographicFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief ConstantDemographicFunction
     *
     *
     */
    class ConstantDemographicFunction : public DemographicFunction {
        
    public:
        ConstantDemographicFunction(const TypedDagNode<double>* t);                                                                                          //!< Default constructor
        ConstantDemographicFunction(const ConstantDemographicFunction &f);                                                                  //!< Copy constructor
        
        // destructor
        virtual                                        ~ConstantDemographicFunction(void);
        ConstantDemographicFunction&                    operator=(const ConstantDemographicFunction &f);                    //!< Assignment operator

        // public methods
        
        // pure virtual public methods
        ConstantDemographicFunction*                    clone(void) const;                                      //!< Clone the ConstantDemographicFunction
        double                                          getDemographic(double t) const;                         //!< demographic function N(t) at time t
        double                                          getIntegral(double start, double finish) const;         //!< Calculates the integral 1/N(x) dx between start and finish.
        
    protected:
        virtual void                                    swapNodeInternal(const DagNode *oldN, const DagNode *newN);

        
    private:
        
        const TypedDagNode<double>*                     theta;

        
    };
    
    // Global functions using the class
    std::ostream&                                       operator<<(std::ostream& o, const ConstantDemographicFunction& x);  //!< Overloaded output stream
    
}

#endif /* ConstantDemographicFunction_H */
