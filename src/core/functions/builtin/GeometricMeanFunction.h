#ifndef GeometricMeanFunction_H
#define GeometricMeanFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Arithmetic mean of a vector of real numbers.
     *
     * This function computes the arithmetic mean of a vector of real number:
     *   y = sum_{i=1}^{n}x[i]/n
     *
     */
    class GeometricMeanFunction : public TypedFunction<double> {
        
    public:
        GeometricMeanFunction(const TypedDagNode< RbVector<double> > * v);
        virtual                                            ~GeometricMeanFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        GeometricMeanFunction*                              clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             vals;
        
    };
    
}


#endif
