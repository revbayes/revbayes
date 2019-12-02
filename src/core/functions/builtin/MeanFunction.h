#ifndef MeanFunction_H
#define MeanFunction_H

#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    
    /**
     * @brief Arithmetic mean of a vector of real numbers.
     *
     * This function computes the arithmetic mean of a vector of real number:
     *   y = sum_{i=1}^{n}x[i]/n
     *
     */
    class MeanFunction : public TypedFunction<double> {
        
    public:
        MeanFunction(const TypedDagNode< RbVector<double> > * v);
        virtual                                            ~MeanFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        MeanFunction*                                       clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             vals;
        
    };
    
}


#endif
