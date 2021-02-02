#ifndef SumFunction_H
#define SumFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Function for computation of the sum of a vector of  real numbers.
     *
     * This class is the function that computes the sum of a vector of numbers.
     * The numbers are passed in as a DAG node whose value type is a std::vector<double>.
     *
     */
    class SumFunction : public TypedFunction<double> {
        
    public:
        SumFunction(const TypedDagNode<RbVector<double> > * v);
        virtual                                            ~SumFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        SumFunction*                                        clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<double> >*              vals;
        
    };
    
}


#endif
