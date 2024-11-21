#ifndef ModuloFunction_H
#define ModuloFunction_H

#include <cstdint>

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
    class ModuloFunction : public TypedFunction<std::int64_t> {
        
    public:
        ModuloFunction(const TypedDagNode<std::int64_t> * l, const TypedDagNode<std::int64_t> *r);
        virtual                                            ~ModuloFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        ModuloFunction*                                     clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                            left;
        const TypedDagNode<std::int64_t>*                            right;
        
    };
    
}


#endif
