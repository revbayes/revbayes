#ifndef StandardDeviationFunction_H
#define StandardDeviationFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief StandardDeviation of a vector of real numbers.
     *
     * This function computes the StandardDeviation of a vector of real number:
     *   y = sum_{i=1}^{n}(x[i]-\bar{x})*(x[i]-\bar{x})/(n-1)
     *
     * This class is the function that computes the standard deviation of a vector of numbers.
     * The numbers are passed in as a DAG node whose value type is a std::vector<double>.
     *
     */

    class StandardDeviationFunction : public TypedFunction<double> {
        
    public:
        StandardDeviationFunction(const TypedDagNode< RbVector<double> > * v);
        virtual                                            ~StandardDeviationFunction(void);                                            //!< Virtual destructor
        
        // public member functions
        StandardDeviationFunction*                          clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             vals;
        
    };
    
}


#endif
