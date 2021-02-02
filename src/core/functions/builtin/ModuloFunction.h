#ifndef ModuloFunction_H
#define ModuloFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Modulo operator function.
     *
     * This function performs the modulo operation, e.g., a % b = c
     * This functions allows this operation to be performed on
     * TypedDagNodes
     *
     */
    class ModuloFunction : public TypedFunction<long> {
        
    public:
        ModuloFunction(const TypedDagNode<long> * l, const TypedDagNode<long> *r);
        virtual                                            ~ModuloFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        ModuloFunction*                                     clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<long>*                            left;
        const TypedDagNode<long>*                            right;
        
    };
    
}


#endif
