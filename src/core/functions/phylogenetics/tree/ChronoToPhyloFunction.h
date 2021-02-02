#ifndef ChronoToPhyloFunction_h
#define ChronoToPhyloFunction_h

#include "Tree.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class ChronoToPhyloFunction : public TypedFunction<Tree> {
        
    public:
        ChronoToPhyloFunction(const TypedDagNode<Tree> *c, const TypedDagNode<RbVector<double> > *r);
        virtual                                            ~ChronoToPhyloFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        ChronoToPhyloFunction*                              clone(void) const;                                                                  //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        // parameters
        const TypedDagNode<Tree>*                           chrono;
        const TypedDagNode<RbVector<double> >*              rates;
    };
    
}

#endif /* ChronoToPhyloFunction_h */
