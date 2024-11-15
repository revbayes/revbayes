//
//  ProductIntegerFunction.hpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 11/15/24.
//

#ifndef ProductIntegerFunction_h
#define ProductIntegerFunction_h


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
    class ProductIntegerFunction : public TypedFunction<long> {
        
    public:
        ProductIntegerFunction(const TypedDagNode<RbVector<long> > * v);
        virtual                      ~ProductIntegerFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        ProductIntegerFunction*      clone(void) const;                                                          //!< Create an independent clone
        void                         update(void);
        
    protected:
        void                         swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<long> >*              vals;
    };
    
}


#endif /* ProductIntegerFunction_h */
