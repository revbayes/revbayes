//
//  ProductFunction.hpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 11/15/24.
//

#ifndef ProductFunction_h
#define ProductFunction_h


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
    class ProductFunction : public TypedFunction<double> {
        
    public:
        ProductFunction(const TypedDagNode<RbVector<double> > * v);
        virtual                                            ~ProductFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        ProductFunction*                                        clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        bool use_safe_product;
        const TypedDagNode<RbVector<double> >*              vals;
    };
    
}


#endif /* ProductFunction_h */
