#ifndef MedianFunction_H
#define MedianFunction_H

#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    
    /**
     * @brief Arithmetic Median of a vector of real numbers.
     *
     * This function computes the arithmetic median of a vector of real number:
     *   y = x[ (x.size()-1)/2 ]                             if x.size() is even
     *   y = (x[ (x.size()-1)/2 ]+x[ x.size()/2 ]) / 2       otherwise
     *
     */
    class MedianFunction : public TypedFunction<double> {
        
    public:
        MedianFunction(const TypedDagNode< RbVector<double> > * v);
        virtual                                            ~MedianFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        MedianFunction*                                     clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             vals;
        
    };
    
}


#endif
