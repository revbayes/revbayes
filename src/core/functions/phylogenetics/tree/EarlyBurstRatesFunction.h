#ifndef EarlyBurstRatesFunction_H
#define EarlyBurstRatesFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class EarlyBurstRatesFunction : public TypedFunction< RbVector<double> > {
        
    public:
        EarlyBurstRatesFunction(const TypedDagNode<Tree> *t, const TypedDagNode<double> *s, const TypedDagNode<double> *l);
        virtual                                            ~EarlyBurstRatesFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        EarlyBurstRatesFunction*                            clone(void) const;                                                                  //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Tree>*                           tau;
        const TypedDagNode<double>*                         sigma;
        const TypedDagNode<double>*                         lambda;
    };
    
}

#endif
